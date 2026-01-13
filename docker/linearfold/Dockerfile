FROM ubuntu:22.04 AS build

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    build-essential \
    python3 \
    ca-certificates \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY external/LinearFold/ /src/
RUN make -j"$(nproc)"

# Runtime: keep python3 for the repository's `linearfold` launcher script.
FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    python3 \
    ca-certificates \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/linearfold
COPY --from=build /src/ /opt/linearfold/

ENV PATH="/opt/linearfold:${PATH}"
ENTRYPOINT ["linearfold"]

