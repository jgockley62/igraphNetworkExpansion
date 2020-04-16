# Igraph Sandbox

### Setup
```{bash}
git clone https://github.com/jgockley62/igraph_Network_Expansion.git

sudo service docker start
sudo usermod -aG docker <USRID>

docker image build -t network ~/igraph_Network_Expansion/Docker/

docker run -v "~/igraph_Network_Expansion/:~/igraph_Network_Expansion/" -e USER=<USERID> -e PASSWORD=<PassWD> -d -p 8787:8787 -p 8080:8080 <ImageID>


```

