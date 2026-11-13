FROM python:3.11.10-slim@sha256:840e180ebcc6e5c8efab209c43f5e40fd2af98cb49db5c7103c90539c56bb30e
WORKDIR /app
COPY . .
ENV PYTHONDONTWRITEBYTECODE=1
# 默认演示保持原有 summary；持续筛查可通过覆盖命令运行，例如：
# docker run --rm adaptive-fdr-monitor python -m analysis_app monitor \
#   --family-weight F0=0.5 --family-weight F1=0.5 --state output/state.json
CMD ["python", "-m", "analysis_app", "summary"]
