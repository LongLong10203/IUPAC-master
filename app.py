from flask import Flask, render_template, request, make_response, redirect, jsonify, session
from main import *
from uuid import uuid4

app = Flask(__name__)
app.secret_key = uuid4().hex

@app.route("/")
def home():
    username = request.cookies.get("username")
    if username is None:
        return redirect("/login")
    return render_template("home.html", username=username)

@app.route("/login", methods=["GET", "POST"])
def login():
    if request.cookies.get("username") is not None:
        return redirect("/")
    if request.method == "POST":
        username = request.form.get("username")
        response = make_response(redirect("/"))
        response.set_cookie("username", username)
        return response
    return render_template("/login.html")

@app.route("/game/getscore")
def get_score():
    return jsonify(score=session["score"])

@app.route("/game/check/<user_answer>")
def check(user_answer: str):
    # print(user_answer, session["answer"])

    if session["answered"]:
        return jsonify(correct=False, cheated=True)

    if user_answer == session["answer"]:
        session["score"] += 1
    session["answered"] = True
    return jsonify(correct=bool(user_answer == session["answer"]))

@app.route("/game")
def game():
    session["answered"] = False
    session["score"] = 0
    session["answer"] = None
    return render_template("/game.html")

@app.route("/game/result")
def game_result():
    return render_template("/result.html")

@app.route("/random_compound")
def random_compound():
    smiles = random_smiles()
    iupac = iupac_name(smiles)
    img_base64 = generate_base64_image(smiles)
    session["answer"] = iupac
    session["answered"] = False
    return jsonify(smiles=smiles, iupac=iupac, img_base64=img_base64)

@app.route("/scoreboard")
def scoreboard():
    return render_template("/scoreboard.html")

@app.route("/logout")
def logout():
    response = make_response(redirect("/login"))
    response.delete_cookie("username")
    return response

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=80, debug=True)