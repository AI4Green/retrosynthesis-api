from flask import Flask

app = Flask(__name__)

from sources.retrosynthesis import retrosynthesis_bp
app.register_blueprint(retrosynthesis_bp)
print("app made")

if __name__ == "__main__":
    app.run()  # to run the app in PyCharm
