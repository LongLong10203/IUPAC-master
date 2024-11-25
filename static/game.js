let score = 0, incorrect = false

function reset() {
    document.getElementById("img").style.display = "none"
    document.getElementById("loader").style.display = "block"
    document.getElementById("answer").value = ""
    document.getElementById("next").style.display = "none"
    document.getElementById("result").innerText = ""
    fetch("/random_compound")
        .then(response => response.json())
        .then(data => {
            document.getElementById("img").src = "data:image/png;base64," + data.img_base64
        })
        .catch(err => {
            console.error(err)
        })
        .then(() => {
            document.getElementById("img").style.display = "block"
            document.getElementById("loader").style.display = "none"
            document.getElementById("answer").disabled = false
        })
}

function update() {
    document.getElementById("score").innerText = score
}

window.addEventListener("DOMContentLoaded", () => {
    reset()
    update()
})

document.getElementById("answer").addEventListener("keydown", (event) => {
    if (event.key === "Enter") {
        fetch(`/game/check/${document.getElementById("answer").value}`)
            .then(res => res.json())
            .then(data => {
                document.getElementById("answer").disabled = true
                document.getElementById("next").style.display = "block"
                if (data.cheated) {
                    document.getElementById("result").innerHTML = "You cheated, right?"
                } else if (data.correct) {
                    document.getElementById("result").innerHTML = "Correct!"
                    ++score
                } else {
                    document.getElementById("result").innerHTML = "Incorrect... The answer is: " + data.answer
                    incorrect = true
                }
                update()
            })
    }
})

function validscore() {
    return fetch("/game/getscore")
        .then(res => res.json())
        .then(data => {
            console.log(data)
            return data.score === score
        })
}

document.getElementById("next").addEventListener("click", () => {
    if (incorrect) {
        validscore().then(valid => {
            if (valid) {
                window.location.href = "/game/result"
            } else {
                alert("You cheated, right?")
                window.location.href = "/"
            }
        })
    } else {
        reset()
    }
})