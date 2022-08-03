# GRAViTy2
A optimised, Python3 implementation of the original GRAViTy (https://github.com/PAiewsakun/GRAViTy)

## Run instructions
1. Install command line requirements (RM < AUTOMATION SCRIPT TO FOLLOW)
1. Navigate to directory: ```cd {path-to-dir}```
1. Install Python libraries with pip: ```pip install -r requirements.txt```
1. *When WP7 finished, initialise Docker container* (RM < TO FOLLOW)
1. *Delete when WP7 finished* Source vars and start API: ```source env_vars.sh && python3 -m uvicorn app.api:app --reload```
1. Navigate to localhost:8000/docs in browser
1. Follow instructions on Swagger UI
