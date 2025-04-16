docker buildx build \
	--platform linux/amd64 \
	--load \
	--no-cache \
	-t us-east1-docker.pkg.dev/wizepair2/wizepair2/wizepair2:v1 \
	--push \
	.
gcloud run deploy wizepair2-v1 \
	--region us-east1 \
	--project wizepair2 \
	--image us-east1-docker.pkg.dev/wizepair2/wizepair2/wizepair2:v1 \
	--platform managed \
	--memory 512Mi \
	--port 56734 \
	--cpu 1 \
	--max-instances 1 \
	--concurrency 1 \
	--timeout 60