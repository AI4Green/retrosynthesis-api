FROM python:3.10-slim

COPY . .
# install the uv package manager
RUN pip install uv
# install app dependencies with uv
RUN uv pip install .


EXPOSE 8000
CMD ["python","src/app.py"]