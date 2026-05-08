 

To build the docker:
go to the file that contraints the Dockerfile and run:

docker build -t arg/cg2o -f Dockerfile ..
docker build --build-arg ENABLE_PARDISO=OFF -t arg/cg2o -f Dockerfile ..
 
to remove the image:
docker rmi  arg/cg2o
 
to remove every thing 
docker system prune -a


To run the container Manually
chmod +x run.sh
./run.sh 
 

After running it 
in vscode 
cntr + alt + p
Dev Containers: Attach to Running Container..
chose the constainer. 

-------------------------------------------------------------------------------------------------------------------------

To run the container using dev container file:
in vscode to add dev container file
Dev Containers: Add Dev Container Configuration Files...
fill the file
Ctrl+Shift+P → "Dev Containers: Reopen in Container"


