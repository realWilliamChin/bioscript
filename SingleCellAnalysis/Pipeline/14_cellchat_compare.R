#!/usr/bin/env Rscript
# Created Time  : 2026/06/29
# Author        : William GoGo
# Pipeline 第 14 步: 把 13 步导出的 *_cellchat.rds 两两 merge, 输出:
#   - 差异交互网络图 (01_*_network.png)
#   - 信号通路秩排序 (02_*_informationFlow.png)
#   - 每个 source/target cluster 的 LR 通讯气泡图 (03_*_dotplot.png)
#   - DEG-based 通讯差异 (04_*_genesExpr_dotplot.png + netMappingDEG.csv)

suppressPackageStartupMessages({
    library(optparse)
    library(CellChat)
    library(patchwork)
    library(dplyr)
    library(stringr)
    library(ggplot2)
})

option_list <- list(
    make_option(c("-i", "--cellchat_dir"), type = "character",
                help = "13 步输出目录, 含 *_cellchat.rds"),
    make_option(c("-o", "--output_dir"), type = "character",
                help = "输出根目录, 自动建 Graph_data/Table_data 子目录"),
    make_option(c("--pairs"), type = "character", default = NULL,
                help = "指定组对, 形如 'HD,Class_III;HD,NoLN'; 默认两两组合"),
    make_option(c("--groups_order"), type = "character", default = NULL,
                help = "可选, 指定 group factor 顺序, 逗号分隔")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$cellchat_dir) || is.null(opt$output_dir)) {
    print_help(OptionParser(option_list = option_list))
    stop("必须提供 --cellchat_dir 与 --output_dir")
}

# ---- 加载所有 cellchat 对象 ----
rds_files <- list.files(opt$cellchat_dir, pattern = "_cellchat\\.rds$", full.names = TRUE)
if (length(rds_files) < 2) stop("至少需要 2 个 *_cellchat.rds 文件")

cellchatLst <- list()
for (f in rds_files) {
    grp <- sub("_cellchat\\.rds$", "", basename(f))
    cellchatLst[[grp]] <- readRDS(f)
    message("加载 ", grp)
}

if (!is.null(opt$groups_order)) {
    order_vec <- trimws(strsplit(opt$groups_order, ",")[[1]])
    cellchatLst <- cellchatLst[order_vec]
}

# ---- 准备组对 ----
if (!is.null(opt$pairs)) {
    pair_list <- lapply(strsplit(opt$pairs, ";")[[1]],
                         function(p) trimws(strsplit(p, ",")[[1]]))
} else {
    all_pairs <- expand.grid(grp1 = names(cellchatLst),
                             grp2 = names(cellchatLst),
                             stringsAsFactors = FALSE)
    all_pairs <- all_pairs[all_pairs$grp1 != all_pairs$grp2, ]
    pair_list <- lapply(split(all_pairs, seq_len(nrow(all_pairs))),
                         function(r) c(r$grp1, r$grp2))
}

dir.create(file.path(opt$output_dir, "Graph_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$output_dir, "Table_data"), recursive = TRUE, showWarnings = FALSE)


# ---- 工具: 单组对的 LR 气泡图 ----
lr_bubble_for_pair <- function(cellchatMergd, grpComp, grp_pair, plotDir, tableDir) {
    sources <- unique(cellchatMergd@idents$joint)
    for (e in sources) {
        targets <- as.character(sources[sources != e])
        for (use_dataset in 1:2) {
            tag <- if (use_dataset == 1) "Incr" else "Decr"
            title_name <- paste0(ifelse(use_dataset == 1, "Increased", "Decreased"),
                                 " signaling in ", grp_pair[1])
            tryCatch({
                gg <- netVisual_bubble(cellchatMergd,
                                       sources.use = e, targets.use = targets,
                                       comparison = c(1, 2), max.dataset = use_dataset,
                                       title.name = title_name,
                                       angle.x = 45, remove.isolate = TRUE,
                                       font.size.title = 14)
                n_labels <- length(unique(gg$data$interaction_name_2))
                height_auto <- max(n_labels * 0.3, 6)
                fname <- paste0("03_source", str_replace_all(e, "/", ""),
                                "_", tag,
                                str_replace_all(grpComp, "_vs_", "_"),
                                "_communicationProbsChange_dotplot.png")
                png(file.path(plotDir, fname),
                    width = 8, height = height_auto, units = "in", res = 300)
                print(gg); dev.off()
            }, error = function(e) message("跳过气泡图: ", e$message))
        }
    }
}


# ---- 工具: 基于 DEG 的 LR 差异通讯 ----
lr_deg_bubble <- function(cellchatMergd, grpComp, grp_pair, plotDir, tableDir) {
    pos.dataset <- grp_pair[1]
    features.name <- paste0(pos.dataset, ".merged")
    tryCatch({
        cellchatMergd <- identifyOverExpressedGenes(cellchatMergd,
                                                    group.dataset = "datasets",
                                                    pos.dataset = pos.dataset,
                                                    features.name = features.name,
                                                    only.pos = FALSE,
                                                    thresh.pc = 0.1, thresh.fc = 0.05)
        net <- netMappingDEG(cellchatMergd, features.name = features.name)
        write.csv(net,
                  file.path(tableDir, paste0(str_replace_all(grpComp, " ", ""),
                                              "_netMappingDEG.csv")),
                  row.names = FALSE, fileEncoding = "UTF-8")

        net_up   <- subsetCommunication(cellchatMergd, net = net,
                                        datasets = grp_pair[1],
                                        ligand.logFC = 0.05, receptor.logFC = NULL)
        net_down <- subsetCommunication(cellchatMergd, net = net,
                                        datasets = grp_pair[1],
                                        ligand.logFC = -0.05, receptor.logFC = NULL)

        sources <- unique(cellchatMergd@idents$joint)
        for (e in sources) {
            targets <- as.character(sources[sources != e])
            for (direction in c("Up", "Down")) {
                lr_use <- if (direction == "Up") net_up else net_down
                if (is.null(lr_use) || nrow(lr_use) == 0) next
                pairLR <- lr_use[, "interaction_name", drop = FALSE]
                title_name <- paste0(direction, "-regulated signaling in ", grp_pair[1])
                tryCatch({
                    gg <- netVisual_bubble(cellchatMergd, pairLR.use = pairLR,
                                           sources.use = e, targets.use = targets,
                                           comparison = c(1, 2),
                                           angle.x = 45, remove.isolate = TRUE,
                                           font.size.title = 14, title.name = title_name)
                    n_labels <- length(unique(gg$data$interaction_name_2))
                    height_auto <- max(n_labels * 0.3, 6)
                    fname <- paste0("04_source", str_replace_all(e, "/", ""),
                                    "_", substr(direction, 1, 4),
                                    str_replace_all(grpComp, "_vs_", "_"),
                                    "_genesExpr_dotplot.png")
                    png(file.path(plotDir, fname),
                        width = 8, height = height_auto, units = "in", res = 300)
                    print(gg); dev.off()
                }, error = function(e) message("跳过 DEG 气泡图: ", e$message))
            }
        }
    }, error = function(e) message("identifyOverExpressedGenes 跳过: ", e$message))
}


# ---- 主循环 ----
for (pair in pair_list) {
    grpComp <- paste0(str_replace_all(pair[1], " ", ""), "_vs_",
                      str_replace_all(pair[2], " ", ""))
    message("==== ", grpComp, " ====")

    objectLst <- cellchatLst[pair]
    cellchatMergd <- mergeCellChat(objectLst, add.names = names(objectLst))
    cellchatMergd@meta$datasets <- factor(cellchatMergd@meta$datasets, levels = pair)

    cells1 <- rownames(cellchatMergd@net[[1]]$weight)
    cells2 <- rownames(cellchatMergd@net[[2]]$weight)
    if (!identical(sort(cells1), sort(cells2))) {
        cellchatMergd <- liftCellChat(cellchatMergd)
    }

    grpCompDir <- file.path(opt$output_dir, "Graph_data", grpComp)
    dir.create(grpCompDir, recursive = TRUE, showWarnings = FALSE)
    tblDir <- file.path(opt$output_dir, "Table_data")

    # 01: diff interaction network
    png(file.path(grpCompDir, paste0("01_", grpComp, "_interactions_network.png")),
        width = 8, height = 6, units = "in", res = 300)
    netVisual_diffInteraction(cellchatMergd, weight.scale = TRUE)
    dev.off()
    png(file.path(grpCompDir, paste0("01_", grpComp, "_strength_network.png")),
        width = 8, height = 6, units = "in", res = 300)
    netVisual_diffInteraction(cellchatMergd, weight.scale = TRUE, measure = "weight")
    dev.off()

    # 02: information flow
    gg1 <- rankNet(cellchatMergd, mode = "comparison", measure = "weight",
                   stacked = TRUE, do.stat = TRUE)
    gg2 <- rankNet(cellchatMergd, mode = "comparison", measure = "weight",
                   stacked = FALSE, do.stat = TRUE)
    ggsave(file.path(grpCompDir, paste0("02_", grpComp, "_informationFlow.png")),
           gg1 + gg2, width = 14, height = 10, dpi = 300)

    # 03 + 04: LR bubble
    lr_bubble_for_pair(cellchatMergd, grpComp, pair, grpCompDir, tblDir)
    lr_deg_bubble(cellchatMergd, grpComp, pair, grpCompDir, tblDir)
}

message("CellChat 比较完成 → ", opt$output_dir)
