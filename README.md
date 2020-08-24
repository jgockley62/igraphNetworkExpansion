# Igraph Sandbox

### Setup
```{bash}
git clone https://github.com/jgockley62/igraph_Network_Expansion.git
cd igraph_Network_Expansion    

#Start Docker
sudo service docker start
sudo usermod -aG docker <USRID>

#Build RStudio Image
docker image build -t network ~/igraph_Network_Expansion/Docker/

#Pull the Cytoscape Network
docker pull biodepot/novnc-cynetworkbma

#Build Caddy image 
cd caddy
docker build -t cyto-caddy .
cd ..

#Hash the password
docker run --rm -it cyto-caddy caddy hash-password -plaintext 'test'

#Change the password to the has in the docker-compose.yaml file"
docker-compose up -d

#RStudio Instance Available at: IP:8787
#Caddy Login to Cytoscale NoVNC VIM Available at: IP:6080


#docker run -v "~/igraph_Network_Expansion/:~/igraph_Network_Expansion/" -e USER=<USERID> -e PASSWORD=<PassWD> -d -p 8787:8787 -p 8080:8080 <ImageID>
#docker pull biodepot/novnc-cynetworkbma
#docker run -v /home/jgockley/igraph_Network_Expansion:/home/jgockley/igraph_Network_Expansion -d -p 80:6080 biodepot/novnc-cynetworkbma
#docker pull cytoscape/ci-igraph
#docker run -v "~/igraph_Network_Expansion/:~/igraph_Network_Expansion/" -e USER=<USERID> -e PASSWORD=<PassWD> -d -p 8080:8080 <IMAGE_ID>

```

