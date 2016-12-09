FROM rocker/shiny

RUN rm /var/lib/apt/lists/http.debian.net_debian_dists_* && \
    apt-get update && \
    apt-get install -y git libxml2-dev libssl-dev ghostscript
 
# Install BST235_project
COPY . /srv/shiny-server/bst235
 
WORKDIR /srv/shiny-server/bst235

# Install packrat packages
RUN Rscript -e 'install.packages("packrat"); \
                packrat::restore()'

# Temporary permissions hack
RUN chown shiny:shiny -R packrat && \
    chown shiny:shiny .gitignore 

# Serve only the scHemeR app (replace /srv/shiny-server with /srv/shiny-server/bst235)
RUN sed -i 's/\/srv\/shiny-server/\/srv\/shiny-server\/bst235/' /etc/shiny-server/shiny-server.conf

# Improve first page load time by not shutting down the R session when idle
RUN sed -i '/location \/ {/a app_idle_timeout 0;' /etc/shiny-server/shiny-server.conf

# Increase app load timeout
RUN sed -i '/location \/ {/a app_init_timeout 60;' /etc/shiny-server/shiny-server.conf

# Add Google Analytics tracking code
# RUN sed -i '/location \/ {/a google_analytics_id UA-37764824-4;' /etc/shiny-server/shiny-server.conf

# Start and expose shiny server
EXPOSE 3838
CMD /usr/bin/shiny-server.sh

