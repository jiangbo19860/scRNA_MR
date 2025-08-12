#!/bin/bash

# 进入项目目录（确保脚本在正确路径执行）
cd /Users/lijiangbo/scRNA_MR/2_src || { echo "目录不存在"; exit 1; }

# 拉取远程最新代码（避免冲突）
git pull origin main

# 暂存所有已跟踪的.R文件（包括修改和删除）
git add -u *.R

# 检查是否有可提交的变化
if git diff --cached --quiet; then
  echo "No changes to commit"
  exit 0
fi

# 提交更改（用日期作为提交信息）
git commit -m "Auto-sync .R scripts: $(date +'%Y-%m-%d %H:%M:%S')"

# 推送到远程仓库
git push origin main
