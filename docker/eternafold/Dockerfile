FROM ubuntu:22.04 AS build

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    build-essential \
    perl \
    ca-certificates \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY external/EternaFold/ /src/
# The repo in `external/` may contain host-built binaries (glibc mismatch). Force a clean rebuild.
RUN make -C /src/src clean || true
RUN make -C /src/src -j"$(nproc)"

FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    ca-certificates \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/eternafold
COPY --from=build /src/parameters /opt/eternafold/parameters
COPY --from=build /src/src/contrafold /usr/local/bin/contrafold
COPY --from=build /src/src/score_prediction /usr/local/bin/score_prediction

ENV ETERNAFOLD_PARAMS="/opt/eternafold/parameters/EternaFoldParams.v1"
ENTRYPOINT ["contrafold"]
