# ----- Builder Stage -----
FROM node:lts AS builder
WORKDIR /app

COPY package*.json ./
RUN npm ci
COPY . .
RUN npm run build

# ----- Runtime Stage -----
FROM node:lts-slim
WORKDIR /app

# Copy only package*.json, then install prod dependencies
COPY package*.json ./
RUN npm ci --production

# Copy over build output and any server files
COPY --from=builder /app/build ./build
COPY --from=builder /app/server.js .

EXPOSE 3000
CMD ["node", "server.js"]
