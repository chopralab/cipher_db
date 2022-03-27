let response = null;
var query;
let plus = document.getElementById("button");
$( "#button" ).click(function() {
plus.innerHTML = '<i class="fa fa-hourglass"></i>';
query = document.getElementById('input-smi').value;
response = null;
add(query);
console.log("adding...");
});

function add(query){
    queryJSON = {
        "smiles": query
    };
    $.ajax({
        url: '/add',
        type: 'POST',
        contentType: 'application/json; charset=utf-8',
        datatype: "json",
        data: JSON.stringify(queryJSON),
        success: function (data) {
            console.log("Recieving search results: " + JSON.stringify(data));
            response = data;
            plus.innerHTML = '<i class="fa fa-check"></i>';
        },
        error: function () {
            console.log('Error');
        }
    });
}

