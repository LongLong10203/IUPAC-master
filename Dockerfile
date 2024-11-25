# Base image
FROM python:3.12

# Set the working directory in the container
WORKDIR /app

# Copy the application files
COPY . .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt
RUN prisma generate

# Expose the port your app runs on
EXPOSE 3016

# Run the Flask application
CMD ["python", "main.py"]