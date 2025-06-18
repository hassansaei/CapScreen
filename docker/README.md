
## Build docker image

```{bash}

docker build --platform linux/amd64  --no-cache -t capscreen:1.0.0.alpha.1 -f docker/DockerFile . 

```

## Run docker image

```{bash}

docker run --platform linux/amd64 -w /opt/capscreen/ -v /Users/hassansaei/Desktop/CapScreen/:/opt/capscreen/input -v /Users/hassansaei/Desktop/CapScreen/:/opt/capscreen/output capscreen:1.0.0.alpha.1 --fastq1 /opt/capscreen/input/P241213-1017ecp-virus_1.fq.gz --fastq2 /opt/capscreen/input/P241213-1017ecp-virus_2.fq.gz --output-dir /opt/capscreen/output --sample-name P241213 --threads 8

```