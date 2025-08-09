# 设置 OpenGWAS令牌: https://api.opengwas.io/profile/,
# https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication
rm(list = ls())
library(ieugwasr)

# 在终端中操作：
# cd /Users/lijiangbo/scRNA_MR
# nano .Renviron
# control + K 删除原来的，然后复制粘贴更新后的，如下面的：
# OPENGWAS_JWT=eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJqaWFuZ2JvMTk4NjBAMTYzLmNvbSIsImlhdCI6MTc1NDU2MDMxMiwiZXhwIjoxNzU1NzY5OTEyfQ.cISh7z7y0nMSaCvyN9rxm42_MCNMc4cHMqIniNpQaen9n4VzLRMIN0MaKbrrxwrd_eDOAQA7mrtev74KUpasEaPuWBMlIKKxV8fWQIUFgdZsHmb-Z9OKopy1_2J0jZH2Jarsm_S-IugvHAgpKTEILzc_BnNtJtAfqI7Sm0pzXJZilC18kYiaR92x3Khuul2ry9Rsgd1IdTDAg-g7XY2TsVungIQspD2xys6kJoqx8-Kvlj-_0C007dos01Y4FLfpeFHTTwSo9EcgxULHNwHbZExZeR515-CfqttIhWBXY_1d-YGL6AxsWo20p1BZZn5CPhWBAXdq-O3_M_W1s12kPQ

# 检查环境变量是否读取成功
Sys.getenv("OPENGWAS_JWT")  # 应返回你的完整令牌字符串

# 检查令牌是否被识别
ieugwasr::get_opengwas_jwt()
ieugwasr::api_status()
ieugwasr::user()

