# 安装/加载 httr 包
install.packages("httr")
library(httr)

# 1. 定义 token 变量，将你的 JWT 令牌粘贴到等号右侧（替换示例内容）
token <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJqaWFuZ2JvMTk4NjBAMTYzLmNvbSIsImlhdCI6MTc1NDU2MDMxMiwiZXhwIjoxNzU1NzY5OTEyfQ.cISh7z7y0nMSaCvyN9rxm42_MCNMc4cHMqIniNpQaen9n4VzLRMIN0MaKbrrxwrd_eDOAQA7mrtev74KUpasEaPuWBMlIKKxV8fWQIUFgdZsHmb-Z9OKopy1_2J0jZH2Jarsm_S-IugvHAgpKTEILzc_BnNtJtAfqI7Sm0pzXJZilC18kYiaR92x3Khuul2ry9Rsgd1IdTDAg-g7XY2TsVungIQspD2xys6kJoqx8-Kvlj-_0C007dos01Y4FLfpeFHTTwSo9EcgxULHNwHbZExZeR515-CfqttIhWBXY_1d-YGL6AxsWo20p1BZZn5CPhWBAXdq-O3_M_W1s12kPQ"  # 这里粘贴你的 token

# 2. 构建请求头（自动拼接 "Bearer " + 令牌）
headers <- add_headers(
  "Authorization" = paste0("Bearer ", token),  # 此处会自动使用上面定义的 token 变量
  "Content-Type" = "application/json"
)
