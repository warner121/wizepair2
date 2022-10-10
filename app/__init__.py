from flask import Flask
#from flask_json_schema import JsonSchema

#from pandarallel import pandarallel
#pandarallel.initialize()

app = Flask(__name__)
#schema = JsonSchema(app)

from app import views
