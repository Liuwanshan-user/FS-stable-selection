# ==========================================
# Stability Selection - 蛋白组特征选择
# ==========================================
# 日期: 2025-10-30
# 版本: 2.0 (优化版)
# 功能: 特征选择 + 核心评估指标
# 改进: 修复参数不一致、增加错误处理、改进代码健壮性
# ==========================================

# ==========================================
# 核心参数设置
# ==========================================
CUTOFF <- 0.65              # 稳定性选择阈值 (修正: 与stabsel保持一致)
B_SAMPLING <- 100           # 抽样次数
FRACTION <- 0.75            # 每次抽样比例
RANDOM_SEED <- 123          # 随机种子
CV_FOLDS <- 5               # 交叉验证折数
PFER <- 2                   # Per-Family Error Rate
INPUT_FILE <- "significant_protein_data.csv"  # 输入文件名 (修正: protein而非metabolite)

# ==========================================
# 安装和加载必要的包
# ==========================================
cat("==========================================\n")
cat("加载R包...\n")
cat("==========================================\n")

packages <- c("stabs", "glmnet", "ggplot2", "dplyr", "pROC", "PRROC", "caret")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("正在安装 %s...\n", pkg))
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}
cat("所有包已加载完成\n\n")

# ==========================================
# 设置工作目录 (改进: 兼容多种运行环境)
# ==========================================
if (require("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  # 在RStudio中运行
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  cat("工作目录: RStudio环境\n")
} else {
  # 在命令行/RScript中运行
  script_dir <- dirname(sys.frame(1)$ofile)
  if (script_dir != "") {
    setwd(script_dir)
  }
  cat("工作目录:", getwd(), "\n")
}

set.seed(RANDOM_SEED)
cat(sprintf("随机种子: %d\n\n", RANDOM_SEED))

# ==========================================
# 1. 读取数据 (增加: 错误处理)
# ==========================================
cat("==========================================\n")
cat("1. 读取数据\n")
cat("==========================================\n")

# 检查文件是否存在
if (!file.exists(INPUT_FILE)) {
  stop(sprintf("错误: 输入文件 '%s' 不存在!\n请确保文件在当前工作目录中: %s",
               INPUT_FILE, getwd()))
}

data <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)
cat("数据维度:", dim(data), "\n")

# 数据结构检查
required_cols <- c("id", "label", "cohort")
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop(sprintf("错误: 数据缺少必需列: %s", paste(missing_cols, collapse = ", ")))
}

# 提取特征和标签
X_raw <- as.matrix(data[, -(1:3)])
y <- as.numeric(data$label)
cohort <- data$cohort
sample_ids <- data$id
feature_names <- colnames(data)[-(1:3)]

cat(sprintf("特征数: %d\n", length(feature_names)))
cat(sprintf("样本数: %d\n", nrow(X_raw)))

# 检查标签是否为二分类
unique_labels <- sort(unique(y))
if (!all(unique_labels %in% c(0, 1))) {
  stop(sprintf("错误: 标签必须是0和1，当前标签: %s",
               paste(unique_labels, collapse = ", ")))
}

# 分离训练集和测试集
train_idx <- which(cohort == "train")
test_idx <- which(cohort == "test")

if (length(train_idx) == 0 || length(test_idx) == 0) {
  stop("错误: 训练集或测试集为空，请检查cohort列的值")
}

# 统计类别分布
cat("\n类别分布:\n")
cat("训练集 - 样本数:", length(train_idx), "\n")
cat("  阴性 (0):", sum(y[train_idx] == 0), "\n")
cat("  阳性 (1):", sum(y[train_idx] == 1), "\n")
cat("测试集 - 样本数:", length(test_idx), "\n")
cat("  阴性 (0):", sum(y[test_idx] == 0), "\n")
cat("  阳性 (1):", sum(y[test_idx] == 1), "\n")

# ==========================================
# 2. 数据预处理 (改进: 增加数据质量检查)
# ==========================================
cat("\n==========================================\n")
cat("2. 数据预处理 (Log2 + 标准化)\n")
cat("==========================================\n")

# 分离训练集和测试集
y_train <- y[train_idx]
y_test <- y[test_idx]

# 数据质量检查
cat("数据质量检查:\n")

# 检查缺失值
missing_count <- sum(is.na(X_raw))
if (missing_count > 0) {
  cat(sprintf("  警告: 发现 %d 个缺失值 (%.2f%%)\n",
              missing_count, missing_count/length(X_raw)*100))
  cat("  处理策略: 用特征中位数填充\n")
  for (i in 1:ncol(X_raw)) {
    na_idx <- is.na(X_raw[, i])
    if (any(na_idx)) {
      X_raw[na_idx, i] <- median(X_raw[!na_idx, i], na.rm = TRUE)
    }
  }
}

# 检查负值和零值
negative_count <- sum(X_raw < 0, na.rm = TRUE)
zero_count <- sum(X_raw == 0, na.rm = TRUE)

cat(sprintf("  负值数量: %d\n", negative_count))
cat(sprintf("  零值数量: %d\n", zero_count))

# Log2转换 (安全处理)
if (negative_count > 0) {
  cat("  警告: 数据包含负值，进行平移处理\n")
  min_val <- min(X_raw, na.rm = TRUE)
  X_raw <- X_raw - min_val + 1  # 平移使所有值 >= 1
}

# 添加小常数避免log(0)
epsilon <- 1e-10
X_log2 <- log2(X_raw + epsilon)

# 检查是否有无限值
if (any(is.infinite(X_log2))) {
  cat("  警告: Log2转换后存在无限值，进行修正\n")
  X_log2[is.infinite(X_log2)] <- NA
  # 用列中位数填充
  for (i in 1:ncol(X_log2)) {
    na_idx <- is.na(X_log2[, i])
    if (any(na_idx)) {
      X_log2[na_idx, i] <- median(X_log2[!na_idx, i], na.rm = TRUE)
    }
  }
}

cat("  数据质量检查完成\n")

# 创建预处理对象（基于训练集）
preProc <- preProcess(
  X_log2[train_idx, ],
  method = c("center", "scale")  # 中心化 + 标准化
)

# 应用预处理（训练集和测试集使用相同参数）
X_train <- predict(preProc, X_log2[train_idx, ])
X_test <- predict(preProc, X_log2[test_idx, ])

cat("预处理完成 (Log2转换 + 中心化 + 标准化)\n")
cat(sprintf("  训练集维度: %d × %d\n", nrow(X_train), ncol(X_train)))
cat(sprintf("  测试集维度: %d × %d\n", nrow(X_test), ncol(X_test)))

# ==========================================
# 3. Stability Selection (修正: 参数一致性)
# ==========================================
cat("\n==========================================\n")
cat("3. Stability Selection\n")
cat("==========================================\n")
cat(sprintf("参数配置:\n"))
cat(sprintf("  Cutoff:   %.2f\n", CUTOFF))
cat(sprintf("  PFER:     %d\n", PFER))
cat(sprintf("  B:        %d\n", B_SAMPLING))
cat(sprintf("  Fraction: %.2f\n", FRACTION))

start_time <- Sys.time()

stab_sel <- stabsel(
  x = X_train,
  y = factor(y_train),
  fitfun = glmnet.lasso,
  cutoff = CUTOFF,              # 修正: 使用定义的参数
  PFER = PFER,                  # 修正: 使用定义的参数
  sampling.type = "SS",
  B = B_SAMPLING,
  assumption = "unimodal",
  args.fitfun = list(fraction = FRACTION)
)

ss_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("完成，用时: %.2f 分钟\n", ss_time))

# 提取稳定性得分
stability_scores <- stab_sel$max
names(stability_scores) <- feature_names

# 选择特征 (修正: 使用一致的阈值)
selected_features <- feature_names[stability_scores >= CUTOFF]
n_selected <- length(selected_features)

cat(sprintf("\n选中特征数: %d / %d (%.1f%%)\n",
            n_selected, length(feature_names),
            n_selected/length(feature_names)*100))

if (n_selected == 0) {
  stop(sprintf("错误: 未选中任何特征，请降低阈值 (当前: %.2f)", CUTOFF))
}

# 显示Top特征
cat("\nTop 10 特征 (按稳定性得分):\n")
top_features <- head(sort(stability_scores, decreasing = TRUE), 10)
for (i in 1:length(top_features)) {
  cat(sprintf("  %2d. %s: %.4f %s\n",
              i, names(top_features)[i], top_features[i],
              ifelse(top_features[i] >= CUTOFF, "[选中]", "")))
}

# ==========================================
# 4. 模型训练和评估 (修正: 使用一致的CV折数)
# ==========================================
cat("\n==========================================\n")
cat("4. 模型训练和评估\n")
cat("==========================================\n")

# 选中特征的数据
X_train_selected <- X_train[, selected_features, drop = FALSE]
X_test_selected <- X_test[, selected_features, drop = FALSE]

# 训练最终模型 (修正: 使用定义的CV_FOLDS)
cat(sprintf("训练Lasso模型 (%d-fold CV)...\n", CV_FOLDS))
final_model <- cv.glmnet(
  X_train_selected,
  y_train,
  family = "binomial",
  alpha = 1,
  type.measure = "auc",
  nfolds = CV_FOLDS              # 修正: 使用定义的参数
)

cat(sprintf("  最优lambda: %.6f\n", final_model$lambda.min))
cat(sprintf("  非零系数数: %d\n",
            sum(coef(final_model, s = "lambda.min") != 0) - 1))  # -1排除截距

# --- 训练集评估 (交叉验证) ---
cat(sprintf("\n训练集评估 (%d-Fold CV):\n", CV_FOLDS))
set.seed(RANDOM_SEED)
folds <- createFolds(y_train, k = CV_FOLDS)

train_cv_pred <- numeric(length(y_train))
for (i in 1:CV_FOLDS) {
  test_idx_cv <- folds[[i]]
  train_idx_cv <- setdiff(1:length(y_train), test_idx_cv)

  cv_model <- cv.glmnet(X_train_selected[train_idx_cv, ], y_train[train_idx_cv],
                        family = "binomial", alpha = 1, type.measure = "auc",
                        nfolds = CV_FOLDS)  # 修正: 使用一致的折数

  train_cv_pred[test_idx_cv] <- as.vector(
    predict(cv_model, newx = X_train_selected[test_idx_cv, , drop = FALSE],
            s = "lambda.min", type = "response")
  )
}

# 训练集指标
train_roc <- roc(y_train, train_cv_pred, quiet = TRUE)
train_auc <- as.numeric(auc(train_roc))
train_pr <- pr.curve(scores.class0 = train_cv_pred[y_train == 1],
                     scores.class1 = train_cv_pred[y_train == 0],
                     curve = TRUE)
train_auprc <- train_pr$auc.integral

# 训练集混淆矩阵
train_pred_class <- ifelse(train_cv_pred > 0.5, 1, 0)
train_cm <- confusionMatrix(factor(train_pred_class, levels = c(0, 1)),
                            factor(y_train, levels = c(0, 1)),
                            positive = "1")

cat(sprintf("  AUC:         %.4f\n", train_auc))
cat(sprintf("  AUPRC:       %.4f\n", train_auprc))
cat(sprintf("  Accuracy:    %.4f\n", train_cm$overall['Accuracy']))
cat(sprintf("  Sensitivity: %.4f\n", train_cm$byClass['Sensitivity']))
cat(sprintf("  Specificity: %.4f\n", train_cm$byClass['Specificity']))

# --- 测试集评估 ---
cat("\n测试集评估:\n")
test_pred <- as.vector(predict(final_model, newx = X_test_selected,
                               s = "lambda.min", type = "response"))

test_roc <- roc(y_test, test_pred, quiet = TRUE)
test_auc <- as.numeric(auc(test_roc))
test_pr <- pr.curve(scores.class0 = test_pred[y_test == 1],
                    scores.class1 = test_pred[y_test == 0],
                    curve = TRUE)
test_auprc <- test_pr$auc.integral

# 测试集混淆矩阵
test_pred_class <- ifelse(test_pred > 0.5, 1, 0)
test_cm <- confusionMatrix(factor(test_pred_class, levels = c(0, 1)),
                           factor(y_test, levels = c(0, 1)),
                           positive = "1")

cat(sprintf("  AUC:         %.4f\n", test_auc))
cat(sprintf("  AUPRC:       %.4f\n", test_auprc))
cat(sprintf("  Accuracy:    %.4f\n", test_cm$overall['Accuracy']))
cat(sprintf("  Sensitivity: %.4f\n", test_cm$byClass['Sensitivity']))
cat(sprintf("  Specificity: %.4f\n", test_cm$byClass['Specificity']))

# 计算过拟合指标
overfit_auc <- train_auc - test_auc
overfit_auprc <- train_auprc - test_auprc

cat("\n过拟合评估:\n")
cat(sprintf("  AUC差异:   %.4f %s\n", overfit_auc,
            ifelse(overfit_auc > 0.1, "[警告: 可能过拟合]", "[良好]")))
cat(sprintf("  AUPRC差异: %.4f %s\n", overfit_auprc,
            ifelse(overfit_auprc > 0.1, "[警告: 可能过拟合]", "[良好]")))

# ==========================================
# 5. 保存结果
# ==========================================
cat("\n==========================================\n")
cat("5. 保存结果\n")
cat("==========================================\n")

if (!dir.exists("results")) {
  dir.create("results")
}

# 5.1 特征稳定性得分
scores_df <- data.frame(
  Feature = feature_names,
  Stability_Score = stability_scores,
  Selected = feature_names %in% selected_features
) %>%
  arrange(desc(Stability_Score))

write.csv(scores_df, "results/feature_stability_scores.csv", row.names = FALSE)
cat("已保存: results/feature_stability_scores.csv\n")

# 5.2 选中的特征列表 (增加: 模型系数)
selected_df <- scores_df %>% filter(Selected == TRUE)

# 提取模型系数
model_coefs <- coef(final_model, s = "lambda.min")
coef_values <- as.vector(model_coefs[-1])  # 排除截距
names(coef_values) <- colnames(X_train_selected)

selected_df$Model_Coefficient <- coef_values[selected_df$Feature]
selected_df <- selected_df %>% arrange(desc(abs(Model_Coefficient)))

write.csv(selected_df, "results/selected_features.csv", row.names = FALSE)
cat(sprintf("已保存: results/selected_features.csv (%d个特征)\n", nrow(selected_df)))

# 5.3 性能指标摘要 (增加: 更多指标)
performance_summary <- data.frame(
  Dataset = c("Train_CV", "Test"),
  N_Samples = c(length(y_train), length(y_test)),
  N_Positive = c(sum(y_train == 1), sum(y_test == 1)),
  N_Negative = c(sum(y_train == 0), sum(y_test == 0)),
  N_Features = c(n_selected, n_selected),
  AUC = c(train_auc, test_auc),
  AUPRC = c(train_auprc, test_auprc),
  Accuracy = c(train_cm$overall['Accuracy'], test_cm$overall['Accuracy']),
  Sensitivity = c(train_cm$byClass['Sensitivity'], test_cm$byClass['Sensitivity']),
  Specificity = c(train_cm$byClass['Specificity'], test_cm$byClass['Specificity']),
  PPV = c(train_cm$byClass['Pos Pred Value'], test_cm$byClass['Pos Pred Value']),
  NPV = c(train_cm$byClass['Neg Pred Value'], test_cm$byClass['Neg Pred Value'])
)

write.csv(performance_summary, "results/performance_summary.csv", row.names = FALSE)
cat("已保存: results/performance_summary.csv\n")

# 5.4 参数配置记录 (新增)
config_df <- data.frame(
  Parameter = c("CUTOFF", "B_SAMPLING", "FRACTION", "CV_FOLDS", "PFER", "RANDOM_SEED"),
  Value = c(CUTOFF, B_SAMPLING, FRACTION, CV_FOLDS, PFER, RANDOM_SEED)
)
write.csv(config_df, "results/analysis_config.csv", row.names = FALSE)
cat("已保存: results/analysis_config.csv\n")

# ==========================================
# 6. 可视化
# ==========================================
cat("\n==========================================\n")
cat("6. 生成可视化\n")
cat("==========================================\n")

# 6.1 ROC曲线
png("results/roc_curves.png", width = 2400, height = 2400, res = 300)
plot(train_roc, col = "steelblue", lwd = 2.5,
     main = "ROC Curves", cex.main = 1.5, cex.lab = 1.2)
plot(test_roc, col = "coral", lwd = 2.5, add = TRUE)
abline(a = 0, b = 1, lty = 2, col = "gray50", lwd = 1.5)
legend("bottomright",
       legend = c(sprintf("Train %d-Fold CV (AUC = %.4f)", CV_FOLDS, train_auc),
                  sprintf("Test (AUC = %.4f)", test_auc)),
       col = c("steelblue", "coral"),
       lwd = 2.5, cex = 1.0, bty = "n")
dev.off()
cat("已保存: results/roc_curves.png\n")

# 6.2 PR曲线
png("results/pr_curves.png", width = 2400, height = 2400, res = 300)
plot(train_pr, col = "steelblue", lwd = 2.5,
     main = "Precision-Recall Curves",
     cex.main = 1.5, cex.lab = 1.2, auc.main = FALSE)
plot(test_pr, col = "coral", lwd = 2.5, add = TRUE)
legend("topright",
       legend = c(sprintf("Train %d-Fold CV (AUPRC = %.4f)", CV_FOLDS, train_auprc),
                  sprintf("Test (AUPRC = %.4f)", test_auprc)),
       col = c("steelblue", "coral"),
       lwd = 2.5, cex = 1.0, bty = "n")
dev.off()
cat("已保存: results/pr_curves.png\n")

# 6.3 混淆矩阵
png("results/confusion_matrices.png", width = 3200, height = 1600, res = 300)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# 训练集混淆矩阵
train_cm_table <- train_cm$table
fourfoldplot(train_cm_table, color = c("coral", "steelblue"),
             conf.level = 0, margin = 1,
             main = sprintf("Train %d-Fold CV\nAcc=%.3f, Sens=%.3f, Spec=%.3f",
                            CV_FOLDS,
                            train_cm$overall['Accuracy'],
                            train_cm$byClass['Sensitivity'],
                            train_cm$byClass['Specificity']))

# 测试集混淆矩阵
test_cm_table <- test_cm$table
fourfoldplot(test_cm_table, color = c("coral", "steelblue"),
             conf.level = 0, margin = 1,
             main = sprintf("Test\nAcc=%.3f, Sens=%.3f, Spec=%.3f",
                            test_cm$overall['Accuracy'],
                            test_cm$byClass['Sensitivity'],
                            test_cm$byClass['Specificity']))

dev.off()
cat("已保存: results/confusion_matrices.png\n")

# 6.4 稳定性得分分布
p_dist <- ggplot(scores_df, aes(x = Stability_Score, fill = Selected)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = CUTOFF, linetype = "dashed",
             color = "red", linewidth = 1.2) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "darkgreen"),
                    labels = c("Not Selected", "Selected")) +
  annotate("text", x = CUTOFF, y = Inf,
           label = sprintf("Cutoff = %.2f", CUTOFF),
           hjust = -0.1, vjust = 2, color = "red",
           size = 5, fontface = "bold") +
  labs(title = "Stability Score Distribution",
       subtitle = sprintf("%d features selected (%.1f%%)",
                          n_selected, n_selected/length(feature_names)*100),
       x = "Stability Score", y = "Count", fill = "") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top")

ggsave("results/stability_score_distribution.png", p_dist,
       width = 10, height = 6, dpi = 300)
cat("已保存: results/stability_score_distribution.png\n")

# 6.5 特征重要性图 (新增)
if (n_selected <= 30) {  # 只显示前30个特征
  top_n <- n_selected
} else {
  top_n <- 30
}

top_features_df <- selected_df %>% head(top_n)

p_importance <- ggplot(top_features_df,
                       aes(x = reorder(Feature, abs(Model_Coefficient)),
                           y = Model_Coefficient)) +
  geom_col(aes(fill = Model_Coefficient > 0), width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "darkgreen", "FALSE" = "coral"),
                    labels = c("Negative", "Positive")) +
  coord_flip() +
  labs(title = sprintf("Top %d Feature Importance", top_n),
       subtitle = "Based on Lasso Model Coefficients",
       x = "", y = "Coefficient", fill = "Direction") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 8))

ggsave("results/feature_importance.png", p_importance,
       width = 10, height = max(6, top_n * 0.25), dpi = 300)
cat("已保存: results/feature_importance.png\n")

# ==========================================
# 7. 完成
# ==========================================
cat("\n==========================================\n")
cat("分析完成!\n")
cat("==========================================\n\n")

cat("核心结果:\n")
cat(sprintf("- 选中特征数: %d / %d (%.1f%%)\n",
            n_selected, length(feature_names),
            n_selected/length(feature_names)*100))
cat(sprintf("- 训练集 AUC: %.4f, AUPRC: %.4f\n", train_auc, train_auprc))
cat(sprintf("- 测试集 AUC: %.4f, AUPRC: %.4f\n", test_auc, test_auprc))
cat(sprintf("- AUC差异: %.4f %s\n", overfit_auc,
            ifelse(overfit_auc > 0.1, "(警告: 可能过拟合)", "(良好)")))
cat(sprintf("- 运行时间: %.2f 分钟\n", ss_time))

cat("\n输出文件:\n")
cat("1. results/feature_stability_scores.csv - 所有特征的稳定性得分\n")
cat("2. results/selected_features.csv - 选中的特征列表+模型系数\n")
cat("3. results/performance_summary.csv - 性能指标汇总\n")
cat("4. results/analysis_config.csv - 分析参数配置 [新增]\n")
cat("5. results/roc_curves.png - ROC曲线\n")
cat("6. results/pr_curves.png - PR曲线\n")
cat("7. results/confusion_matrices.png - 混淆矩阵\n")
cat("8. results/stability_score_distribution.png - 稳定性得分分布\n")
cat("9. results/feature_importance.png - 特征重要性图 [新增]\n")
cat("==========================================\n")
