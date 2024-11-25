#!/bin/sh

prisma db push

# This will exec the CMD from your Dockerfile
exec "$@"