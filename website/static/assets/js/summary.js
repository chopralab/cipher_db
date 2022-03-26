var resultsList = document.getElementById('results');
var response = null;
let query = window.location.href;
console.log(query);
query = query.split('/')[4];
console.log(query);
search(query);
console.log("searching...");


$( "#copy-json" ).click(function() {
    navigator.clipboard.writeText(JSON.stringify(response));
    icon = document.getElementById('clip-check');
    icon.classList.remove('fa-clipboard');
    icon.classList.add('fa-check');
});

function capitalizeFirstLetter(string) {
  return string.charAt(0).toUpperCase() + string.slice(1);
}

function search(query){
    queryJSON = {
        "inchikey": query
    };
    $.ajax({
        url: '/info',
        type: 'POST',
        contentType: 'application/json; charset=utf-8',
        datatype: "json",
        data: JSON.stringify(queryJSON),
        success: function (data) {
            console.log("Recieving search results: " + JSON.stringify(data));
            response = data;
            console.log(response)
            renderResults(response);
            //document.getElementById("loading").remove();
        },
        error: function () {
            console.log('Error');
        }
    });
}

function generateData(cols,bounds){
    let vals = [];
    for(let i = 0; i < cols; i++){
        vals.push(Math.random() * (bounds.max - bounds.min) + bounds.min);
    }
    return vals;
}
var options = {
      series: [
      {
        name: 'Naloxone',
        data: generateData(9, {
          min: -30,
          max: 55
        })
      },
      {
        name: 'Desired biosignature',
        data: generateData(9, {
          min: -30,
          max: 55
        })
      }
    ],
      chart: {
      height: 350,
      type: 'heatmap',
    },
    plotOptions: {
      heatmap: {
        shadeIntensity: 0.5,
        radius: 0,
        useFillColorAsStroke: true,
        colorScale: {
          ranges: [{
              from: -30,
              to: 5,
              name: 'Partial agonism',
              color: '#00A100'
            },
            {
              from: 6,
              to: 20,
              name: 'Agonism',
              color: '#128FD9'
            },
            {
              from: 21,
              to: 45,
              name: 'Antagonism',
              color: '#FFB200'
            },
            {
              from: 46,
              to: 55,
              name: 'No effect',
              color: '#050505'
            }
          ]
        }
      }
    },
    dataLabels: {
      enabled: false
    },
    stroke: {
      width: 1
    },
    title: {
      text: ''
    },
     xaxis: {
         type: 'category',
         tickPlacement: 'between',
         axisTicks: {
            show: true
         },
  tickAmount: 9,
  range: 9,
  labels: {
    show: true,
    rotate: -90,
    rotateAlways: true,
    formatter: function(value,timestamp,opts){
        console.log(value);
        return ["MOR","DOR","KOR","NOR","D2DD","D2LDR","D3DR","NMDAR","AMPAR"][value-1]
    }
  }
}
    };

    var chart = new ApexCharts(document.querySelector("#sig"), options);
    chart.render();


function renderResults(results){
    let svg = document.getElementById('cmpd-svg');
    svg.innerHTML = results.props.svg;
    let pubchem = document.getElementById('pubchem');
    let rdkit = document.getElementById('rdkit');
    phtml = '<p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">PubChem<br></p>';
    phtml += '<table style="table-layout: fixed; color: grey;border-style: solid;" class="table"><thead><tr><th scope="col">Property Name</th><th scope="col">Property Value</th></tr></thead><tbody>';
    rhtml = '<p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">RDKit<br></p>';
    rhtml += '<table style="table-layout: fixed; color: grey;border-style: solid;" class="table"><thead><tr><th scope="col">Property Name</th><th scope="col">Property Value</th></tr></thead><tbody>';
    for(var key of Object.keys(results.props.pubchem)){
        if (key == 'modified'){
            continue;
        }
        phtml += '<tr><td>'+key+'</td><td style="word-wrap:break-word;">'+results.props.pubchem[key]+'</td></tr>';
    }
    phtml += '</tbody></table>';
    pubchem.innerHTML = phtml;
    for(var key of Object.keys(results.props.rdkit)){
        if (key == 'modified'){
            continue;
        }
        rhtml += '<tr><td>'+key+'</td><td style="word-wrap:break-word;">'+results.props.rdkit[key]+'</td></tr>';
    }
    rhtml += '</tbody></table>';
    rdkit.innerHTML = rhtml;
    
    let assaysContent = document.getElementById('assays-content');
    ahtml = "";
    console.log(results.assays)
    for(let i=0; i<results.assays.length; i++){
        ahtml += '<p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Assay  '+results.assays[i]._id+'<br></p>';
        ahtml += '<table style="table-layout: fixed; color: grey;border-style: solid;" class="table"><thead><tr><th scope="col">Property Name</th><th scope="col">Property Value</th></tr></thead><tbody>';
        for (var key of Object.keys(results.assays[i])){
            if (key == "_id"){
                continue;
            }
            ahtml += '<tr><td>'+key+'</td><td style="word-wrap:break-word;">'+results.assays[i][key]+'</td></tr>';
        }
        ahtml += '</tbody></table>';
    }
    assaysContent.innerHTML = ahtml;
}
