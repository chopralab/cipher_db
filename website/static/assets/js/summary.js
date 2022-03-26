var resultsList = document.getElementById('results');
var response = null;
let query = window.location.href;
console.log(query);
query = query.split('/')[4];
console.log(query);
search(query);
console.log("searching...");


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
    rhtml = '<p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">RDKit<br></p>';
    for(var key of Object.keys(results.props.pubchem)){
        if (key == 'modified'){
            continue;
        }
        phtml += '<p data-bs-toggle="tooltip" data-bss-tooltip="" data-bs-placement="left" style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+key+': '+results.props.pubchem[key]+'<br></p>';
    }
    pubchem.innerHTML = phtml;
    for(var key of Object.keys(results.props.rdkit)){
        if (key == 'modified'){
            continue;
        }
        rhtml += '<p data-bs-toggle="tooltip" data-bss-tooltip="" data-bs-placement="left" style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+key+': '+results.props.rdkit[key]+'<br></p>';
    }
    rdkit.innerHTML = rhtml;
}
