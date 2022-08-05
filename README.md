```
   _____ _____       __      ___ _______     __      _____  
  / ____|  __ \     /\ \    / (_)__   __|    \ \    / /__ \ 
 | |  __| |__) |   /  \ \  / / _   | |_   _   \ \  / /   ) |
 | | |_ |  _  /   / /\ \ \/ / | |  | | | | |   \ \/ /   / / 
 | |__| | | \ \  / ____ \  /  | |  | | |_| |    \  /   / /_ 
  \_____|_|  \_\/_/    \_\/   |_|  |_|\__, |     \/   |____|
                                       __/ |                
                                      |___/                 
 ```                                                                         
GRAViTy V2 adaptation by Mayne, R., Aiewsakun, P., Simmonds., P. et al. (2022)

Based on original GRAViTy software https://github.com/PAiewsakun/GRAViTy

Aiewsakun, P., Simmonds, P. The genomic underpinnings of eukaryotic virus taxonomy: creating a sequence-based framework for family-level virus classification. Microbiome 6, 38 (2018). https://doi.org/10.1186/s40168-018-0422-7

## Run instructions
1. Install command line requirements (RM < AUTOMATION SCRIPT TO FOLLOW)
1. Navigate to directory: ```cd {path-to-dir}```
1. Install Python libraries with pip: ```pip install -r requirements.txt```
1. *When WP7 finished, initialise Docker container* (RM < TO FOLLOW)
1. *Delete when WP7 finished* Source vars and start API: ```source env_vars.sh && python3 -m uvicorn app.api:app --reload```
1. Navigate to localhost:8000/docs in browser
1. Follow instructions on Swagger UI
