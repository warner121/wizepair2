# ---- Stage 1: Builder ----
FROM python:3.12-slim AS builder

# Install system dependencies and Poetry
ENV POETRY_VERSION=1.7.1 \
    POETRY_HOME=/opt/poetry \
    PATH="/opt/poetry/bin:$PATH"

RUN apt-get update && \
    apt-get install -y --no-install-recommends curl build-essential && \
    rm -rf /var/lib/apt/lists/* && \
    curl -sSL https://install.python-poetry.org | python3 -

WORKDIR /app

# Copy only what's needed to install dependencies
COPY pyproject.toml poetry.lock ./
RUN poetry config virtualenvs.create false && \
    poetry install --only main --no-interaction --no-ansi

# ---- Stage 2: Runtime ----
FROM python:3.12-slim

# Set environment early to persist Poetry PATH
ENV POETRY_HOME=/opt/poetry
ENV PATH="${POETRY_HOME}/bin:${PATH}"

WORKDIR /app

# Copy Poetry from builder
COPY --from=builder /opt/poetry /opt/poetry

# Copy dependencies and app code
COPY --from=builder /usr/local/lib/python3.12/site-packages /usr/local/lib/python3.12/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin
COPY . .

# Set environment variables
ENV PORT=8000 \
    WORKERS=1 \
    TIMEOUT=60 \
    HOST=0.0.0.0 \
    PYTHONUNBUFFERED=1

# Run as non-root user
RUN useradd -m appuser && chown -R appuser /app
USER appuser

EXPOSE ${PORT}

# Entrypoint using Uvicorn for FastAPI
CMD ["sh", "-c", "exec uvicorn main:app --host=${HOST} --port=${PORT} --workers=${WORKERS} --timeout-keep-alive=${TIMEOUT}"]