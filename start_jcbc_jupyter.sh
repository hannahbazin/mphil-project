#!/bin/bash

REMOTE_USER="hmlb2"
REMOTE_HOST="mti-ai-srv-01.jcbc.private.cam.ac.uk"
REMOTE_PORT=8088
REMOTE_DIR="/home/hmlb2/mphil-project"
VENV_PATH="bat_proximity_analysis/bin/activate"

# Run command in interactive login shell via -t and -l
ssh -t ${REMOTE_USER}@${REMOTE_HOST} "bash -l -c 'cd ${REMOTE_DIR} && source ${VENV_PATH} && jupyter notebook --no-browser --port=${REMOTE_PORT}'"

