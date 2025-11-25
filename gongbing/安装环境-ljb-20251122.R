#https://www.yuque.com/post_gwas/qtlmr/szamugaqzgznhd6h
# 查看简要版本信息
R.version.string
setwd("D:\\1_Ljb\\QTLMR\\r_env\\QTLMR\\QTLMR")
install.packages("renv")
renv::restore()



#如果使用的是R语言4.3.x及4.4.x版本 在百度网盘 “R环境配置”文件夹下有个“QTLMR.zip"文件，
#已经预装好所有依赖包，解压后可使用下面函数激活
#独立R环境下使用QTLMR包，解决安装依赖包的麻烦，其它版本R语言则不支持！


##安装QTLMRget辅助安装包##
install.packages("remotes")
remotes::install_github("Hortoner/QTLMRget")

# QTLMRget::Activate_QTLMR(QTLMR_path = "D:\\1_Ljb\\QTLMR")


##在线安装或更新QTLMR包##
QTLMRget::install_QTLMR()

#安装Rtools

# install_Rtools()                   #或者在已下载的"软件"文件夹中手动安装

#植入IEU OPENGWAS token,需要自行申请，有效期15天

install_token(OPENGWAS_JWT =" ")      #设置完成后需要重启R

ieugwasr::user()                      #验证IEU OPEN GWAS token有效期

usethis::edit_r_environ()             #查看全局环境设置


setwd("D:/1_Ljb/gongbing")
library(QTLMR)

#添加conda目录权限，所有用户拥有所有权限才可以执行，没有权限可能会导入告警
python_envs_install(envs_path = "D:\\ProgramData\\miniconda3\\envs",
                    python_path = "D:\\1_Ljb\\QTLMR\\qtlmr_python_env")
