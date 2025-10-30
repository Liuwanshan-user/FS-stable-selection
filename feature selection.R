# ==========================================
# Stability Selection - 蛋白组/代谢组特征选择
# ==========================================
# 日期: 2025-10-30
# 功能: 特征选择 + 核心评估指标
# ==========================================

# ==========================================
# 核心参数设置
# ==========================================
STABSEL_CUTOFF <- 0.65      # stabsel函数的cutoff参数 (必须 > 0.5)
FINAL_CUTOFF <- 0.50        # 最终筛选特征的阈值 (可以 < 0.5)
B_SAMPLING <- 100           # 抽样次数
FRACTION <- 0.75            # 每次抽样比例
RANDOM_SEED <- 123          # 随机种子
CV_FOLDS <- 5               # 交叉验证折数
PFER <- 2                   # Per-Family Error Rate
INPUT_FILE <- "significant_metabolite_data.csv"  # 输入文件名

# ==========================================
# 安装和加载必要的包
# ==========================================
cat("加载R包...\n")
packages <- c("caret", "stabs", "glmnet", "ggplot2", "dplyr", "pROC", "PRROC")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(RANDOM_SEED)

# ==========================================
# 1. 读取数据
# ==========================================
cat("\n==========================================\n")
cat("1. 读取数据\n")
cat("==========================================\n")

data <- read.csv(INPUT_FILE)
cat("数据维度:", dim(data), "\n")

# 提取特征和标签
X_raw <- as.matrix(data[, -(1:3)])
y <- as.numeric(data$label)
cohort <- data$cohort
sample_ids <- data$id
feature_names <- colnames(data)[-(1:3)]

# 分离训练集和测试集
train_idx <- which(cohort == "train")
test_idx <- which(cohort == "test")

# 统计类别分布
cat("\n类别分布:\n")
cat("训练集 - 样本数:", length(train_idx), "\n")
cat("  阴性 (0):", sum(y[train_idx] == 0), "\n")
cat("  阳性 (1):", sum(y[train_idx] == 1), "\n")
cat("测试集 - 样本数:", length(test_idx), "\n")
cat("  阴性 (0):", sum(y[test_idx] == 0), "\n")
cat("  阳性 (1):", sum(y[test_idx] == 1), "\n")

# ==========================================
# 2. 数据预处理
# ==========================================
cat("\n==========================================\n")
cat("2. 数据预处理 (Log2 + 标准化)\n")
cat("==========================================\n")

# 分离训练集和测试集
y_train <- y[train_idx]
y_test <- y[test_idx]

# Log2转换
X_log2 <- log2(X_raw)

# 创建预处理对象（基于训练集）
preProc <- preProcess(
  X_log2[train_idx, ],
  method = c("center", "scale")  # 中心化 + 标准化
)

# 应用预处理（训练集和测试集使用相同参数）
X_train <- predict(preProc, X_log2[train_idx, ])
X_test <- predict(preProc, X_log2[test_idx, ])

cat("预处理完成 (使用caret::preProcess)\n")

# ==========================================
# 3. Stability Selection
# ==========================================
cat("\n==========================================\n")
cat("3. Stability Selection\n")
cat("==========================================\n")
cat(sprintf("参数: stabsel_cutoff=%.2f, final_cutoff=%.2f, B=%d, fraction=%.2f\n",
            STABSEL_CUTOFF, FINAL_CUTOFF, B_SAMPLING, FRACTION))

start_time <- Sys.time()

stab_sel <- stabsel(
  x = X_train,
  y = factor(y_train),
  fitfun = glmnet.lasso,
  cutoff = STABSEL_CUTOFF,     # stabsel函数的cutoff (必须 > 0.5)
  PFER = PFER,
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

# 选择特征（使用FINAL_CUTOFF，可以 < 0.5）
selected_features <- feature_names[stability_scores >= FINAL_CUTOFF]
n_selected <- length(selected_features)

cat(sprintf("\n选中特征数: %d / %d (%.1f%%)\n",
            n_selected, length(feature_names),
            n_selected/length(feature_names)*100))

if (n_selected == 0) {
  stop("错误: 未选中任何特征，请降低FINAL_CUTOFF阈值")
}

# ==========================================
# 4. 模型训练和评估
# ==========================================
cat("\n==========================================\n")
cat("4. 模型训练和评估\n")
cat("==========================================\n")

# 选中特征的数据
X_train_selected <- X_train[, selected_features, drop = FALSE]
X_test_selected <- X_test[, selected_features, drop = FALSE]

# 训练最终模型
cat("训练Lasso模型...\n")
final_model <- cv.glmnet(
  X_train_selected,
  y_train,
  family = "binomial",
  alpha = 1,
  type.measure = "auc",
  nfolds = CV_FOLDS
)

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
                        nfolds = CV_FOLDS)

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

cat(sprintf("  AUC: %.4f\n", train_auc))
cat(sprintf("  AUPRC: %.4f\n", train_auprc))
cat(sprintf("  Accuracy: %.4f\n", train_cm$overall['Accuracy']))
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

cat(sprintf("  AUC: %.4f\n", test_auc))
cat(sprintf("  AUPRC: %.4f\n", test_auprc))
cat(sprintf("  Accuracy: %.4f\n", test_cm$overall['Accuracy']))
cat(sprintf("  Sensitivity: %.4f\n", test_cm$byClass['Sensitivity']))
cat(sprintf("  Specificity: %.4f\n", test_cm$byClass['Specificity']))

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

# 5.2 选中的特征列表
selected_df <- scores_df %>% filter(Selected == TRUE)
write.csv(selected_df, "results/selected_features.csv", row.names = FALSE)
cat(sprintf("已保存: results/selected_features.csv (%d个特征)\n", nrow(selected_df)))

# 5.3 性能指标摘要
performance_summary <- data.frame(
  Dataset = c("Train_CV", "Test"),
  N_Samples = c(length(y_train), length(y_test)),
  N_Positive = c(sum(y_train == 1), sum(y_test == 1)),
  N_Negative = c(sum(y_train == 0), sum(y_test == 0)),
  AUC = c(train_auc, test_auc),
  AUPRC = c(train_auprc, test_auprc),
  Accuracy = c(train_cm$overall['Accuracy'], test_cm$overall['Accuracy']),
  Sensitivity = c(train_cm$byClass['Sensitivity'], test_cm$byClass['Sensitivity']),
  Specificity = c(train_cm$byClass['Specificity'], test_cm$byClass['Specificity'])
)

write.csv(performance_summary, "results/performance_summary.csv", row.names = FALSE)
cat("已保存: results/performance_summary.csv\n")

# 5.4 导出选中特征的表达数据
selected_data <- data.frame(
  id = sample_ids,
  label = y,
  cohort = cohort
)

# 添加选中特征的原始表达数据
selected_data <- cbind(selected_data, X_raw[, selected_features, drop = FALSE])

write.csv(selected_data, "results/selected_features_expression_data.csv", row.names = FALSE)
cat(sprintf("已保存: results/selected_features_expression_data.csv (id + label + cohort + %d个特征)\n", n_selected))

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
  geom_vline(xintercept = FINAL_CUTOFF, linetype = "dashed",
             color = "red", linewidth = 1.2) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "darkgreen"),
                    labels = c("Not Selected", "Selected")) +
  annotate("text", x = FINAL_CUTOFF, y = Inf,
           label = sprintf("Cutoff = %.2f", FINAL_CUTOFF),
           hjust = -0.1, vjust = 2, color = "red",
           size = 5, fontface = "bold") +
  labs(title = "Stability Score Distribution",
       subtitle = sprintf("%d features selected", n_selected),
       x = "Stability Score", y = "Count", fill = "") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top")

ggsave("results/stability_score_distribution.png", p_dist,
       width = 10, height = 6, dpi = 300)
cat("已保存: results/stability_score_distribution.png\n")

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
cat(sprintf("- 运行时间: %.1f 分钟\n", ss_time))

cat("\n输出文件:\n")
cat("1. results/feature_stability_scores.csv - 所有特征的稳定性得分\n")
cat("2. results/selected_features.csv - 选中的特征列表\n")
cat("3. results/performance_summary.csv - 性能指标汇总\n")
cat("4. results/selected_features_expression_data.csv - 选中特征的表达数据 (id+label+cohort+特征)\n")
cat("5. results/roc_curves.png - ROC曲线\n")
cat("6. results/pr_curves.png - PR曲线\n")
cat("7. results/confusion_matrices.png - 混淆矩阵\n")
cat("8. results/stability_score_distribution.png - 稳定性得分分布\n")
cat("==========================================\n")
