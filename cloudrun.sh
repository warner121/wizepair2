 docker build -t gcr.io/wizepair2/wizepair2:v1 -f Dockerfile .
 docker push gcr.io/wizepair2/wizepair2:v1
 gcloud run deploy wizepair2-v1 --image gcr.io/wizepair2/wizepair2:v1 --region us-east1 --platform managed --memory 512Mi --port 56734 --cpu 2 --max-instances 15 --timeout 600 --concurrency 6
