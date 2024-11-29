from flask import Flask, render_template, request, make_response, redirect, jsonify, session
from ext_chem import *
from uuid import uuid4
import bcrypt
from prisma import Prisma
from datetime import datetime

app = Flask(__name__)
app.secret_key = uuid4().hex
# disable cache
app.config["SEND_FILE_MAX_AGE_DEFAULT"] = 0

prisma = Prisma()

async def connect():
    if not prisma.is_connected():
        await prisma.connect()

async def disconnect():
    if prisma.is_connected():
        await prisma.disconnect()

@app.route("/")
def home():
    if not session["logged_in"]:
        return redirect("/login")
    
    return render_template("home.html", username=session["username"])

@app.route("/login", methods=["GET", "POST"])
def login():
    if session["logged_in"]:
        return redirect("/")
    
    if request.method == "POST":
        username = request.form.get("username")
        session["username"] = username
        return redirect("/login/password")
    
    return render_template("/login.html")

@app.route("/login/password", methods=["GET", "POST"])
async def login_password():
    if session["logged_in"]:
        return redirect("/")
    
    username = session["username"]
    if username is None: # user directly jumps to this page instead of going through the /login page
        return redirect("/login")
    
    await connect()
    user = await prisma.user.find_unique(
        where={
            "username": username
        }
    )
    await disconnect()

    if request.method == "POST":
        password = request.form.get("password")

        if user is None:
            await connect()
            await prisma.user.create(
                data={
                    "username": username,
                    "password": bcrypt.hashpw(password.encode("utf-8"), bcrypt.gensalt()).decode("utf-8")
                }
            )
            await disconnect()
        else: # user already exists
            # check password
            if not bcrypt.checkpw(password.encode("utf-8"), user.password.encode("utf-8")):
                return render_template("/login_password.html", new=False, wrong_password=True)

        session["logged_in"] = True
        return redirect("/")
    
    return render_template("/login_password.html", new=(user is None))

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
    return jsonify(correct=bool(user_answer == session["answer"]), cheated=False, answer=session["answer"])

@app.route("/game")
def game():
    session["answered"] = False
    session["score"] = 0
    session["answer"] = None
    return render_template("/game.html")

@app.route("/game/result")
async def game_result():
    username = session["username"]

    await connect()
    user = await prisma.user.find_unique(
        where={"username": username}
    )
    new_high_score = session["score"] > user.max_score
    if new_high_score:
        await connect()
        await prisma.user.update(
            where={
                "username": username
            },
            data={
                "max_score": session["score"],
                "updated_at": datetime.now()
            }
        )
    await disconnect()

    return render_template("/result.html", new_high_score=new_high_score)

@app.route("/random_compound")
def random_compound():
    smiles = random_smiles()
    iupac = iupac_name(smiles)
    img_base64 = generate_base64_image(smiles)
    session["answer"] = iupac
    session["answered"] = False
    return jsonify(img_base64=img_base64)

@app.route("/scoreboard")
async def scoreboard():
    await connect()
    users = await prisma.user.find_many(
        order={"max_score": "desc"}
    )
    await disconnect()
    return render_template("/scoreboard.html", users=users, current_user=session["username"])

@app.route("/logout")
def logout():
    session["logged_in"] = False
    session["username"] = None
    return redirect("/login")

@app.route("/account/delete")
async def delete_account():
    await connect()
    await prisma.user.delete(
        where={"username": session["username"]}
    )
    await disconnect()
    return make_response(redirect("/logout"))

@app.before_request
async def before_request():
    if session.get("logged_in") is None:
        session["logged_in"] = False
    if session.get("username") is None:
        session["username"] = None

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=3016, debug=True)