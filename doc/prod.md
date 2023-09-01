Docker up

~~~
cd container
docker compose --project-name gwas2eqtl_pleiotropy_prod --env-file env_prod -f docker-compose.yml up --build --force-recreate --remove-orphans -d
~~~

For the following like dev.md, ....
