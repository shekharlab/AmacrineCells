# packages
require(xgboost)
source("xgboost_script.R")
xg_data = readRDS("XgDataMacaqueBC.rds")

# Get matrices and cell ids. Note that the matrices represent log(TPM+1) gene expression values
train_data = xg_data$train_mat
test_data = xg_data$test_mat
train_id = xg_data$train_id
test_id = xg_data$test_id

# We have pre-selected variable genes for each dataset using the method described in Pandey, Shekhar et al., Current Biology, 2018. Here, we use the common set of variable genes for training
var.genes_train = intersect(xg_data$train_vargenes, xg_data$test_vargenes)

# subset matrices to only include variable genes as features
train_data = train_data[var.genes_train,]; test_data = test_data[var.genes_train,]

# Train xgboost model
bst_model = XGBoost_train(train_data, train_labels = train_id, var.genes = var.genes_train, do.scale = TRUE)

# Use trained model to predict on test data
numberOfClasses = length(levels(train_id))
test_data_scale = xgb.DMatrix(scale(t(test_data), scale=bst_model$scale_var, center = bst_model$scale_mean ))
test_pred <- predict(bst_model$bst_model, newdata = test_data_scale)
test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                          ncol=length(test_pred)/numberOfClasses)

# Find best class for each cell
test_pred_margins = apply(test_prediction,2,max)
test_predlabels=apply(test_prediction,2,which.max)
names(test_predlabels) = colnames(test_data)
test_pred_names = levels(train_id)[test_predlabels]
names(test_pred_names) = names(test_predlabels)

# Plot confusion matrix
C = table(test_id, test_pred_names)
plotConfusionMatrix(C, order="Row", xlab.use = "Training clusters", ylab.use = "Test clusters", title.use = "Performance Confusion Matrix")
