FROM rocker/tidyverse:latest

RUN apt-get update\
&& apt-get install -y libpython-dev\
&& apt-get install -y libpython3-dev\
&& apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev libglpk-dev\
&& apt-get install -y curl nano less libcurl4-openssl-dev gawk\
&& apt-get install -y build-essential chrpath libssl-dev libxft-dev\
&& apt-get install -y libfreetype6 libfreetype6-dev\
&& apt-get install -y libfontconfig1 libfontconfig1-dev\
&& export PHANTOM_JS="phantomjs-1.9.8-linux-x86_64"\
&& wget https://bitbucket.org/ariya/phantomjs/downloads/$PHANTOM_JS.tar.bz2\
&& tar xvjf $PHANTOM_JS.tar.bz2\
&& mv $PHANTOM_JS /usr/local/share\
&& ln -sf /usr/local/share/$PHANTOM_JS/bin/phantomjs /usr/local/bin\
&& apt-get install -y php php-common gcc\
&& apt-get install -y imagemagick\
&& apt-get install -y php-imagick\
&& apt-get install -y git\
&& apt-get install -y tabix\
&& apt-get install -y bcftools\
&& apt-get install -y awscli\
&& apt-get install -y less\
&& apt-get install -y vsftpd\
&& apt-get install -y vcftools\
&& apt-get install -y coreutils\
&& apt install -y python-pip\
&& pip install synapseclient\
&& pip install aws-mfa\
&& pip install awsmfa\
&& git clone https://github.com/jgockley62/igraph_Network_Expansion.git \
&& R -e "install.packages('BiocManager')"\
&& R -e "BiocManager::install('biomaRt')"\ 
#&& R -e "devtools::install_github('th1vairam/CovariateAnalysis@dev')"\ 
&& R -e "devtools::install_github('brian-bot/githubr')"\ 
&& R -e "devtools::install_github('blogsdon/spike/spike')"\ 
&& R -e "install.packages('doParallel')"\
&& R -e "install.packages('data.table')"\
&& R -e "install.packages('dplyr')"\
&& R -e "install.packages('foreach')"\
&& R -e "install.packages('ggplot2')"\
&& R -e "install.packages('bcv')"\
&& R -e "install.packages('knitr')"\
&& R -e "install.packages('plyr')"\
&& R -e "install.packages('rmarkdown')"\
&& R -e "install.packages('reticulate')"\
&& R -e "install.packages('plotrix', repos = 'http://cran.us.r-project.org')"\
&& R -e "install.packages('reticulate', repos = 'http://cran.us.r-project.org')"\
&& R -e "install.packages('optmatch', repos = 'http://cran.us.r-project.org')" \
&& R -e "install.packages('pheatmap', repos = 'http://cran.us.r-project.org')"\
&& R -e "BiocManager::install('scater')"\
&& R -e "install.packages('tidyr')"\
&& R -e "install.packages('glmnet', repos = 'http://cran.us.r-project.org')"\
&& R -e "BiocManager::install('sva')"\
&& R -e "BiocManager::install('WGCNA')"\
&& R -e "BiocManager::install('cqn')"\
&& R -e "devtools::install_github('rlbarter/superheat')"\
&& R -e "devtools::install_github('blogsdon/spike/spike')"\
&& R -e "install.packages('ggrepel', repos = 'http://cran.us.r-project.org')"\
&& R -e "install.packages('kableExtra', repos = 'http://cran.us.r-project.org')"\
&& R -e "BiocManager::install('limma')"\
&& R -e "install.packages('vbsr', repos='http://cran.us.r-project.org')"\
&& R -e "BiocManager::install('edgeR')"\
&& R -e "BiocManager::install('ComplexHeatmap')"


RUN apt-get install -y r-cran-rjava \
&& R -e "install.packages('igraph')" \
&& R -e "install.packages('reshape')" \
&& R -e "install.packages('RColorBrewer')"\
&& R -e "install.packages('parallel')" \
&& R -e "BiocManager::install('DOSE', ask=F)"\
&& R -e "BiocManager::install('org.Hs.eg.db')"

RUN sudo apt-get install -y software-properties-common \
&& sudo apt-add-repository -y ppa:marutter/rrutter \
#&& sudo apt-get -y update \
&& sudo apt-get -y upgrade \
&& sudo apt-get -y install g++ \
&& sudo apt-get -y install libcurl4-openssl-dev

RUN  sudo apt-get -y install liblzma-dev \
&& sudo apt-get -y install libbz2-dev \
&& sudo apt-get -y install libpcre++-dev \
&& sudo apt-get -y install openjdk-11-jdk \ 
&& sudo apt-get -y install libxml2-dev \
&& sudo R CMD javareconf \
&& R -e "BiocManager::install('paxtoolsr')" \
&& R -e "BiocManager::install('clusterProfiler')" \
&& R -e "install.packages('RJSONIO')"\
&& R -e "BiocManager::install('RCy3')"\
&& sudo apt-get install -y glpk-utils libglpk-dev glpk-doc \
&& R -e "install.packages('Rglpk')" \
&& R -e "install.packages('corpcor')" \
&& R -e  "devtools::install_github('jhrcook/HotNetvieweR')"
