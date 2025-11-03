# 设置项目主目录
project_dir <- "~/scRNAseq"

# 定义子文件夹路径（按层级结构）
subfolders <- c(
  "data",
  "data/raw",
  "data/processed",
  "data/metadata",
  "scripts",
  "scripts/utils",
  "results",
  "results/figures",
  "results/tables",
  "results/reports",
  "docs",
  "config",
  "notebooks"
)

# 批量创建文件夹
for (folder in subfolders) {
  full_path <- file.path(project_dir, folder)
  if (!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE)
    cat(paste0("创建文件夹: ", full_path, "\n"))
  } else {
    cat(paste0("文件夹已存在: ", full_path, "\n"))
  }
}

# 创建.gitignore文件（空文件）
gitignore_path <- file.path(project_dir, ".gitignore")
if (!file.exists(gitignore_path)) {
  file.create(gitignore_path)
  cat("创建.gitignore文件\n")
}

# 创建README.md文件（基础说明）
readme_path <- file.path(project_dir, "docs", "README.md")
if (!file.exists(readme_path)) {
  readme_content <- "### scMR_COVID19_Analysis 项目\n\n"
  readme_content <- paste(readme_content, "单细胞孟德尔随机化分析COVID-19的免疫机制与药物靶点\n\n", sep = "\n")
  readme_content <- paste(readme_content, "## 项目结构\n", sep = "\n")
  readme_content <- paste(readme_content, "```\nscMR_COVID19_Analysis/\n├── data/              # 原始数据与处理后数据\n", sep = "\n")
  readme_content <- paste(readme_content, "├── scripts/           # R脚本\n", sep = "\n")
  readme_content <- paste(readme_content, "├── results/           # 分析结果\n", sep = "\n")
  readme_content <- paste(readme_content, "├── docs/              # 文档\n", sep = "\n")
  readme_content <- paste(readme_content, "├── config/            # 配置文件\n", sep = "\n")
  readme_content <- paste(readme_content, "└── .gitignore         # 版本控制忽略文件\n```", sep = "\n")

  writeLines(readme_content, readme_path)
  cat("创建README.md文件\n")
}

cat("项目文件夹构建完成！\n")
