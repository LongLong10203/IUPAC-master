from flask import Flask, render_template, request, make_response, redirect, jsonify, session
from ext_chem import *
from uuid import uuid4
import bcrypt
from prisma import Prisma

app = Flask(__name__)
app.secret_key = uuid4().hex

prisma = Prisma()

@app.route("/")
def home():
    if request.cookies.get("logged_in") is None:
        return redirect("/login")
    return render_template("home.html", username=request.cookies.get("username"))

@app.route("/login", methods=["GET", "POST"])
def login():
    if request.cookies.get("logged_in") is not None:
        return redirect("/")
    
    if request.method == "POST":
        username = request.form.get("username")
        response = make_response(redirect("/login/password"))
        response.set_cookie("username", username)
        return response
    
    return render_template("/login.html")

@app.route("/login/password", methods=["GET", "POST"])
async def login_password():
    if request.cookies.get("logged_in"):
        return redirect("/")
    
    username = request.cookies.get("username")
    if username is None: # user directly jumps to this page instead of going through the /login page
        return redirect("/login")
    
    await prisma.connect()
    user = await prisma.user.find_first(
        where={
            "username": username
        }
    )
    await prisma.disconnect()

    if request.method == "POST":
        password = request.form.get("password")

        if user is None:
            await prisma.connect()
            await prisma.user.create(
                data={
                    "username": username,
                    "password": bcrypt.hashpw(password.encode("utf-8"), bcrypt.gensalt()).decode("utf-8")
                }
            )
            await prisma.disconnect()
        else: # user already exists
            # check password
            if not bcrypt.checkpw(password.encode("utf-8"), user.password.encode("utf-8")):
                return render_template("/login_password.html", new=False, wrong_password=True)

        response = make_response(redirect("/"))
        response.set_cookie("logged_in", "true")
        return response
    
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
    username = request.cookies.get("username")

    await prisma.connect()
    user = await prisma.user.find_first(
        where={"username": username}
    )
    new_high_score = session["score"] > user.max_score
    if new_high_score:
        await prisma.user.update(
            where={"username": username},
            data={"max_score": session["score"]}
        )
    await prisma.disconnect()

    return render_template("/result.html", new_high_score=new_high_score)

@app.route("/random_compound")
def random_compound():
    smiles = random_smiles()
    iupac = iupac_name(smiles)
    img_base64 = generate_base64_image(smiles)
    session["answer"] = iupac
    session["answered"] = False
    return jsonify(smiles=smiles, iupac=iupac, img_base64=img_base64)

@app.route("/scoreboard")
async def scoreboard():
    await prisma.connect()
    users = await prisma.user.find_many(
        order={"max_score": "desc"}
    )
    await prisma.disconnect()
    return render_template("/scoreboard.html", users=users, current_user=request.cookies.get("username"))

@app.route("/logout")
def logout():
    response = make_response(redirect("/login"))
    response.delete_cookie("username")
    response.delete_cookie("logged_in")
    return response

@app.route("/account/delete")
async def delete_account():
    await prisma.connect()
    await prisma.user.delete(
        where={"username": request.cookies.get("username")}
    )
    await prisma.disconnect()
    return make_response(redirect("/logout"))

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=3016, debug=True)