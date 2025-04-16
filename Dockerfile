# ---- Stage 1: Builder ----
FROM python:3.11-slim AS builder

# Install system dependencies and Poetry
ENV POETRY_VERSION=1.7.1 \
    POETRY_HOME=/opt/poetry \
    PATH="/opt/poetry/bin:$PATH"

RUN apt-get update && \
    apt-get install -y --no-install-recommends curl build-essential && \
    rm -rf /var/lib/apt/lists/* && \
    curl -sSL https://install.python-poetry.org | python3 -

WORKDIR /app

# Copy only what’s needed to install dependencies
COPY pyproject.toml poetry.lock ./
RUN poetry config virtualenvs.create false && \
    poetry install --only main --no-interaction --no-ansi

# ---- Stage 2: Runtime ----
FROM python:3.11-slim

WORKDIR /app

# Copy dependencies and binaries from builder
COPY --from=builder /usr/local/lib/python3.11/site-packages /usr/local/lib/python3.11/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin

# Copy app code
COPY . .

# Set environment variables
ENV PORT=5000 \
    WORKERS=1 \
    THREADS=1 \
    WORKER_CLASS=gthread \
    TIMEOUT=60 \
    PYTHONUNBUFFERED=1

# Run as non-root user
RUN useradd -m appuser && chown -R appuser /app
USER appuser

EXPOSE ${PORT}

# Entrypoint
CMD ["sh", "-c", "exec gunicorn main:app --bind=0.0.0.0:${PORT} --workers=${WORKERS} --threads=${THREADS} --worker-class=${WORKER_CLASS} --timeout=${TIMEOUT}"]
