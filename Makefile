ifneq (,$(wildcard ./.deploy.env))
        include ./.deploy.env
endif

# ==================================================================================== #
# Deploy
# ==================================================================================== #
## deploy/staging: deploy web to staging environment
.PHONY: deploy
deploy: put-docker-compose build
	@ssh -i $(PROD_SSH_KEY_PATH) -p $(PROD_PORT) -t $(PROD_USER)@$(PROD_HOST) " \
		cd $(PROD_DIR) && \
		docker compose -f docker-compose.yaml pull && \
		docker compose -f docker-compose.yaml up -d --no-deps --remove-orphans && \
		docker system prune -a --volumes --force"

.PHONY: build
build:
	DOCKER_DEFAULT_PLATFORM=linux/amd64 docker compose build && \
	docker compose push

.PHONY: put-docker-compose
deploy/put-docker-compose:
	@echo "uploading docker-compose.yaml"
	@ssh -i $(PROD_SSH_KEY_PATH) -p $(PROD_PORT) $(PROD_USER)@$(PROD_HOST) "mkdir -p $(PROD_DIR)"
	@scp -i $(PROD_SSH_KEY_PATH) -P $(PROD_PORT) ./docker-compose.yaml $(PROD_USER)@$(PROD_HOST):$(PROD_DIR)/docker-compose.yaml
	@echo "docker-compose.yaml uploaded"
