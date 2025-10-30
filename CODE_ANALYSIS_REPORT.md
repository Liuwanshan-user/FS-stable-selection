# 代码审查与优化报告

## 项目信息
- **项目名称**: Stability Selection 蛋白组特征选择
- **原始文件**: feature selection.R
- **优化文件**: feature_selection_optimized.R
- **审查日期**: 2025-10-30
- **版本**: 2.0

---

## 一、发现的问题总结

### 🔴 严重错误 (Critical Issues)

#### 1. 参数不一致导致逻辑错误
**位置**: 第11行、第106行、第122行

**原始代码**:
```r
# 第11行: 定义参数
CUTOFF <- 0.50

# 第106行: stabsel函数中硬编码不同的值
stab_sel <- stabsel(
  cutoff = 0.65,  # ❌ 未使用定义的CUTOFF参数
  ...
)

# 第122行: 使用原定义的参数筛选
selected_features <- feature_names[stability_scores >= CUTOFF]  # 使用0.50
```

**问题分析**:
- Stability Selection在stabsel中使用0.65作为cutoff
- 但在选择特征时使用0.50作为阈值
- 导致选择逻辑不一致，可能选中不应选中的特征

**修复方案**:
```r
CUTOFF <- 0.65  # 统一阈值

stab_sel <- stabsel(
  cutoff = CUTOFF,  # ✅ 使用定义的参数
  ...
)

selected_features <- feature_names[stability_scores >= CUTOFF]  # 一致性
```

---

#### 2. 数据预处理缺少安全检查
**位置**: 第78行

**原始代码**:
```r
X_log2 <- log2(X_raw)  # ❌ 直接转换，未检查负值/零值
```

**问题分析**:
- 如果数据包含负值，`log2(负值)` = NaN
- 如果数据包含零值，`log2(0)` = -Inf
- 导致后续分析全部失败

**修复方案**:
```r
# 检查负值和零值
negative_count <- sum(X_raw < 0, na.rm = TRUE)
zero_count <- sum(X_raw == 0, na.rm = TRUE)

if (negative_count > 0) {
  min_val <- min(X_raw, na.rm = TRUE)
  X_raw <- X_raw - min_val + 1  # 平移使所有值 >= 1
}

# 添加小常数避免log(0)
epsilon <- 1e-10
X_log2 <- log2(X_raw + epsilon)

# 检查无限值
if (any(is.infinite(X_log2))) {
  # 用列中位数填充
  ...
}
```

---

#### 3. 环境依赖问题
**位置**: 第30行

**原始代码**:
```r
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # ❌ 只能在RStudio运行
```

**问题分析**:
- 只能在RStudio中运行
- 命令行执行 `Rscript feature_selection.R` 会报错
- 限制了代码的可移植性

**修复方案**:
```r
if (require("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  # 在RStudio中运行
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  # 在命令行/RScript中运行
  script_dir <- dirname(sys.frame(1)$ofile)
  if (script_dir != "") {
    setwd(script_dir)
  }
}
```

---

### 🟡 次要问题 (Minor Issues)

#### 4. 交叉验证折数不一致
**位置**: 第15行、第152行

**原始代码**:
```r
# 第15行
CV_FOLDS <- 5

# 第152行
final_model <- cv.glmnet(..., nfolds = 10)  # ❌ 硬编码为10
```

**修复**: 统一使用 `nfolds = CV_FOLDS`

---

#### 5. 文件名与功能不匹配
**位置**: 第40行

**原始代码**:
```r
data <- read.csv("significant_metabolite_data.csv")  # ❌ metabolite(代谢物)
```

**问题**: 注释说"蛋白组特征选择"，但文件名是metabolite(代谢物)

**修复**: 改为 `"significant_protein_data.csv"` 或使用可配置参数

---

#### 6. 缺少文件存在性检查
**原始代码**: 直接读取CSV，如果文件不存在会报错

**修复**:
```r
if (!file.exists(INPUT_FILE)) {
  stop(sprintf("错误: 输入文件 '%s' 不存在!", INPUT_FILE))
}
```

---

#### 7. 缺少缺失值检查
**问题**: 未检查数据是否有NA值

**修复**:
```r
missing_count <- sum(is.na(X_raw))
if (missing_count > 0) {
  cat(sprintf("警告: 发现 %d 个缺失值\n", missing_count))
  # 用中位数填充
  for (i in 1:ncol(X_raw)) {
    na_idx <- is.na(X_raw[, i])
    if (any(na_idx)) {
      X_raw[na_idx, i] <- median(X_raw[!na_idx, i], na.rm = TRUE)
    }
  }
}
```

---

#### 8. 重复加载包
**位置**: 第21行、第75行

```r
# 第21行: 加载caret
packages <- c("stabs", "glmnet", "ggplot2", "dplyr", "pROC", "PRROC", "caret")

# 第75行: 重复加载
library(caret)  # ❌ 已经在第21-28行加载过
```

**修复**: 删除第75行的重复加载

---

## 二、代码改进内容

### ✅ 改进1: 参数管理优化
- 所有参数在开头集中定义
- 全局一致使用，避免硬编码
- 新增 `INPUT_FILE` 参数方便配置

### ✅ 改进2: 增强数据质量检查
- 缺失值检测与处理
- 负值/零值检测与处理
- 无限值检测与修正
- 数据结构验证（必需列检查）

### ✅ 改进3: 改善错误处理
- 文件存在性检查
- 标签格式验证（必须是0/1）
- 训练集/测试集非空检查
- 未选中特征时的明确错误提示

### ✅ 改进4: 增强可移植性
- 兼容RStudio和命令行环境
- 工作目录自动检测
- 更好的环境适配性

### ✅ 改进5: 新增功能
1. **过拟合评估**
   ```r
   overfit_auc <- train_auc - test_auc
   cat(sprintf("AUC差异: %.4f %s\n", overfit_auc,
               ifelse(overfit_auc > 0.1, "[警告: 可能过拟合]", "[良好]")))
   ```

2. **模型系数输出**
   - 在选中特征CSV中添加模型系数列
   - 按系数绝对值排序

3. **特征重要性可视化**
   - 新增 `feature_importance.png`
   - 展示Top特征的Lasso系数

4. **参数配置记录**
   - 新增 `analysis_config.csv`
   - 记录所有分析参数，便于复现

5. **更详细的输出信息**
   - Top 10特征显示
   - 数据质量报告
   - 模型非零系数统计

### ✅ 改进6: 代码可读性
- 增加详细注释说明修改位置
- 统一代码风格
- 清晰的分隔符和输出格式

---

## 三、性能指标对比

### 原始代码输出:
- 特征稳定性得分
- 性能指标(AUC, AUPRC等)
- 4个可视化图表

### 优化代码输出:
- 所有原始输出 ✅
- **新增**: 模型系数
- **新增**: 过拟合评估
- **新增**: 参数配置记录
- **新增**: 特征重要性图
- **新增**: 数据质量报告

---

## 四、使用建议

### 参数调优建议:

1. **CUTOFF (稳定性阈值)**
   - 默认: 0.65
   - 范围: 0.5-0.9
   - 越高越保守，特征越少

2. **PFER (Per-Family Error Rate)**
   - 默认: 2
   - 范围: 1-5
   - 控制假阳性率

3. **B_SAMPLING (抽样次数)**
   - 默认: 100
   - 推荐: 100-500
   - 更多抽样更稳定但更慢

4. **FRACTION (抽样比例)**
   - 默认: 0.75
   - 范围: 0.5-0.9
   - Stability Selection标准做法

5. **CV_FOLDS (交叉验证折数)**
   - 默认: 5
   - 推荐: 5或10
   - 小数据集用5，大数据集用10

---

## 五、测试清单

使用优化代码前，请确保:

- [ ] 输入文件名正确（修改 `INPUT_FILE` 参数）
- [ ] 数据格式正确（必需列: id, label, cohort）
- [ ] 标签为0/1二分类
- [ ] 有足够的训练集和测试集样本
- [ ] R版本 >= 4.0.0
- [ ] 所有依赖包已安装

---

## 六、关键修复对照表

| 问题类型 | 原始代码 | 优化代码 | 影响 |
|---------|---------|---------|------|
| 参数不一致 | cutoff=0.65 vs CUTOFF=0.50 | 统一使用CUTOFF=0.65 | 🔴 严重 |
| Log2转换 | 未检查负值/零值 | 完整的安全检查 | 🔴 严重 |
| 环境依赖 | 仅RStudio | 兼容命令行 | 🔴 严重 |
| CV折数 | 硬编码10 | 使用CV_FOLDS | 🟡 次要 |
| 文件检查 | 无 | 增加存在性检查 | 🟡 次要 |
| 缺失值处理 | 无 | 中位数填充 | 🟡 次要 |

---

## 七、文件对比

### 原始文件: `feature selection.R`
- 行数: 369
- 输出文件: 7个
- 错误处理: 基本

### 优化文件: `feature_selection_optimized.R`
- 行数: 556 (+187行)
- 输出文件: 9个 (+2个)
- 错误处理: 完善
- 新增功能: 5项

---

## 八、总结

### 主要改进:
1. ✅ 修复了3个严重逻辑错误
2. ✅ 修复了5个次要问题
3. ✅ 增加了完善的数据质量检查
4. ✅ 增强了代码健壮性和可移植性
5. ✅ 新增了5个实用功能
6. ✅ 保持了原有代码的核心逻辑

### 建议:
- **立即使用优化版本**，避免参数不一致导致的错误结果
- **保留原始版本**作为备份
- **根据实际数据调整参数**（特别是CUTOFF和PFER）
- **检查过拟合指标**，如果AUC差异>0.1需要重新调参

---

**审查人员**: Claude Code
**优化完成时间**: 2025-10-30
