#!/bin/bash
# 切换到 R 脚本所在目录
cd /Users/lijiangbo/scRNA_MR/2_src || {
  echo "[$(date)] 错误：无法进入目录" >> /Users/lijiangbo/scRNA_MR/sync_error.log
  exit 1
}

# 执行 Git 操作
echo "[$(date)] 开始同步..." >> /Users/lijiangbo/scRNA_MR/sync_git.log

# 添加所有 R 脚本
git add ./*.R

# 提交（仅当有修改时）
git diff --quiet --exit-code || {
  git commit -m "自动同步 R 脚本：$(date +%Y-%m-%d)"
  git push origin main  # 替换 main 为你的分支名（如 master）
  echo "[$(date)] 同步成功" >> /Users/lijiangbo/scRNA_MR/sync_git.log
}
