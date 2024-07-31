from flask import Blueprint

retrosynthesis_bp = Blueprint('retrosynthesis', __name__)

from sources.retrosynthesis import routes
