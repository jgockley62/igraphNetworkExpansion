# Igraph Sandbox

### Install Git, Docker, and Docker-Compose if needed
This is configured for AWS EC-2 instance
```{bash}
sudo yum install -y git
git version

sudo amazon-linux-extras install docker
docker version

sudo curl -L https://github.com/docker/compose/releases/download/1.22.0/docker-compose-$(uname -s)-$(uname -m) -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose
docker-compose version

```

### Setup The Git Repo
```{bash}
git clone https://github.com/jgockley62/igraph_Network_Expansion.git
cd igraph_Network_Expansion    
```

### Start Docker
```{bash}
sudo service docker start
sudo usermod -aG docker <USR_ID>
```

### Build RStudio Images
```{bash}
#RStudio
docker image build -t network ~/igraph_Network_Expansion/Docker/

#Cytoscape Linux VIM
docker pull biodepot/novnc-cynetworkbma

#Caddy proxy https login 
docker build -t cyto-caddy caddy/.
```
### Build the envronment object containing login credentials
```{bash}
#Hash the password
docker run --rm -it cyto-caddy caddy hash-password -plaintext 'test'
```

### Build Containers
```{bash}

docker-compose up -d

```

### Ported Browser Access
RStudio Instance Available at: https://<AWS Instance IP>:8787

Caddy Proxy Login to Access Cytoscape NoVNC VIM Available at: https://<AWS Instance IP>:6080

### Shut Containers Down
```{bash}

docker-compose down -v

```

### Deprecated code
```{bash}

docker run -v "~/igraph_Network_Expansion/:~/igraph_Network_Expansion/" -e USER=<USERID> -e PASSWORD=<PassWD> -d -p 8787:8787 <ImageID>

docker run -v /home/jgockley/igraph_Network_Expansion:/home/jgockley/igraph_Network_Expansion -d -p 6080:6080 biodepot/novnc-cynetworkbma

```

