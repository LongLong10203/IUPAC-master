document.getElementById("logout-btn").addEventListener("click", () => {
    if (confirm("Are you sure you want to logout?")) {
        window.location.href = "/logout"
    }
})

document.getElementById("acc-del-btn").addEventListener("click", () => {
    let username = document.cookie.split(';').filter(item => item.trim().startsWith('username=')).reduce((prev, current) => {
        return current.split('=')[1]
    }, '')
    if (confirm(`Are you sure you want to delete your account (${username}), including losing the max score?\nThis process is irreversible.`)) {
        window.location.href = "/account/delete"
    }
})