from flask import Flask
from flask_cors import CORS

from execute_models import execute_models_api

app = Flask(__name__)
app.register_blueprint(execute_models_api)
CORS(app)




@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"

@app.route('/aircadia')
def aaa():
    return "Hello, AirCADia!"



if __name__ == '__main__':
    app.run(debug=True, port=3005)