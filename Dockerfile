# Multi-stage build for optimized image size
# Stage 1: Builder
FROM python:3.11-slim AS builder

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    curl \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install Poetry
ENV POETRY_VERSION=1.7.1
ENV POETRY_HOME=/opt/poetry
RUN curl -sSL https://install.python-poetry.org | python3 -
ENV PATH="${POETRY_HOME}/bin:${PATH}"

WORKDIR /app

# Copy dependency specifications
COPY pyproject.toml poetry.lock ./

# Install project dependencies
RUN poetry config virtualenvs.create false && \
    poetry install --only main --no-interaction --no-ansi

# Stage 2: Runtime
FROM python:3.11-slim

WORKDIR /app

# Copy installed dependencies from builder
COPY --from=builder /usr/local/lib/python3.11/site-packages /usr/local/lib/python3.11/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin

# Copy application code
COPY . .

# Environment variables
ENV PORT=5000
ENV PYTHONUNBUFFERED=1
ENV GUNICORN_CMD_ARGS="--bind=0.0.0.0:${PORT} --workers=1 --threads=1 --worker-class=gthread --timeout 60"

# Security: Run as non-root user
RUN useradd -m appuser && chown -R appuser /app
USER appuser

EXPOSE ${PORT}

# Start command
CMD ["gunicorn", "main:app"]
