To build the docker, run : 
> chmod +x run_pipeline.sh
> ./run_pipeline.sh build config.yaml

or the classic : 
> docker build . blablabla



To run the docker, run :
> ./run_pipeline.sh run config.yaml

To save the docker image (if you later want to export it) : 
> docker image save -o predimel-deconv.tar predimel-deconv:latest

And then to load the .tar to another machine's docker :
> docker image load -i predimel-deconv.tar
