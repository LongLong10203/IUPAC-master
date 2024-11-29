let score = 0, incorrect = false, answered = false

function reset() {
    document.getElementById("img").style.display = "none"
    document.getElementById("loader").style.display = "block"
    document.getElementById("answer").value = ""
    document.getElementById("submit").style.display = "block"
    document.getElementById("save-and-exit").style.display = "block"
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
            document.getElementById("answer").focus()
            answered = false
        })
}

function update() {
    document.getElementById("score").innerText = score
}

window.addEventListener("DOMContentLoaded", () => {
    reset()
    update()
})

function check_ans() {
    fetch(`/game/check/${document.getElementById("answer").value.trim().toLowerCase()}`)
        .then(res => res.json())
        .then(data => {
            document.getElementById("answer").disabled = true
            document.getElementById("submit").style.display = "none"
            document.getElementById("save-and-exit").style.display = "none"
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
            answered = true
        })
}

document.getElementById("answer").addEventListener("keydown", (event) => {
    if (event.key === "Enter") {
        check_ans()
    }
})

document.getElementById("submit").addEventListener("click", () => {
    check_ans()
})

function validscore() {
    return fetch("/game/getscore")
        .then(res => res.json())
        .then(data => {
            // console.log(data)
            return data.score === score
        })
}

function next() {
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
}

document.getElementById("save-and-exit").addEventListener("click", () => {
    if (confirm("Are you sure you want to save and exit?")) {
        incorrect = true
        next()
    }
})

document.getElementById("next").addEventListener("click", () => {
    next()
})

document.addEventListener("keydown", (event) => {
    if (event.key === "Enter" && answered) {
        event.preventDefault()
        next()
    }
})