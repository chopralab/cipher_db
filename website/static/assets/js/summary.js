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

function renderBiosig(cmpdName,inchikey,cmpdBiosig,desired){
    bioList = [];
    desList = [];
    bioList.push(cmpdBiosig[1].AMPAR);
    bioList.push(cmpdBiosig[1].D2LDR);
    bioList.push(cmpdBiosig[1].DRD2);
    bioList.push(cmpdBiosig[1].DRD3);
    bioList.push(cmpdBiosig[1].NMDAR);
    bioList.push(cmpdBiosig[1].deltaOR);
    bioList.push(cmpdBiosig[1].kappaOR);
    bioList.push(cmpdBiosig[1].muOR);
    bioList.push(cmpdBiosig[1].nociceptinOR);
    desList.push(desired[1].AMPAR);
    desList.push(desired[1].D2LDR);
    desList.push(desired[1].DRD2);
    desList.push(desired[1].DRD3);
    desList.push(desired[1].NMDAR);
    desList.push(desired[1].deltaOR);
    desList.push(desired[1].kappaOR);
    desList.push(desired[1].muOR);
    desList.push(desired[1].nociceptinOR);
        
    var options = {
      series: [
      {
        name: cmpdName,
        data: bioList
      },
      {
        name: 'Desired biosignature',
        data: desList
      }
    ],
      chart: {
      height: 350,
      type: 'heatmap',
    },
    plotOptions: {
      heatmap: {
        shadeIntensity: 1,
        radius: 0,
        useFillColorAsStroke: true,
        colorScale: {
          ranges: [{
              from: 0,
              to: 1,
              name: 'Partial agonism',
              color: '#00A100'
            },
            {
              from: 1,
              to: 2,
              name: 'Agonism',
              color: '#128FD9'
            },
            {
              from: 2,
              to: 3,
              name: 'Antagonism',
              color: '#FFB200'
            },
            {
              from: 3,
              to: 4,
              name: 'Unknown effect',
              color: '#050505'
            },
            {
              from: -9999999999999999999,
              to: -1,
              name: 'Null',
              color: '#ed3e3e'
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
      text: capitalizeFirstLetter(cmpdName)
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
        //console.log(value);
        return ["AMPAR","D2LDR","DRD2","DRD3","NMDAR","deltaOR","kappaOR","muOR","nociceptinOR"][value-1]
        
    }
  }
}
    };

    var chart = new ApexCharts(document.querySelector("#sig"), options);
    chart.render();
}



function renderResults(results){
    let cmpdName = results.ids[0].name;
    console.log(cmpdName);
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
    if (results.assays.length == 0) {
        ahtml = '<p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">No assay information available<br></p>';
    }
    assaysContent.innerHTML = ahtml;
    
    let synth = document.getElementById('synth-routes');
    let synthsList = "";
    //for (let j=0; j<results.synths.length; j++){
    for (let j=0; j<Math.min(5,results.synths.length); j++){
        synthsList += '<p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Retrosynthetic route '+String(j+1)+'<br></p><img src="data:image/png;base64,'+results.synths[j]+'"><hr>'
    }
    if (results.synths.length == 0) {
        synthsList = '<p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">No synthetic information available<br></p>';
    }
    synth.innerHTML = synthsList;
    
    let biosig = results.biosigs;
    let desired = results.desired;
    console.log(biosig);
    console.log(desired);
    for (var key of Object.keys(desired[0])){
        if (desired[1][key] == -1){
            continue;
        }
        if (desired[0][key] == "partial agonist"){
            desired[1][key] = desired[1][key];
        }
        if (desired[0][key] == "binder"){
            desired[1][key] = desired[1][key] + 3;
        }
        if (desired[0][key] == "agonism"){
            desired[1][key] = desired[1][key] + 1;
        }
        if (desired[0][key] == "antagonism"){
            desired[1][key] = desired[1][key] + 2;
        }
        if (desired[0][key] == "unknown effect"){
            desired[1][key] = desired[1][key] + 3;
        }
    }
    for (var key of Object.keys(biosig[0])){
        if (desired[1][key] == -1){
            continue;
        }
        if (biosig[0][key] == "partial agonist"){
            biosig[1][key] = biosig[1][key];
        }
        if (biosig[0][key] == "binder"){
            biosig[1][key] = biosig[1][key] + 3;
        }
        if (biosig[0][key] == "agonism"){
            biosig[1][key] = biosig[1][key] + 1;
        }
        if (biosig[0][key] == "antagonism"){
            biosig[1][key] = biosig[1][key] + 2;
        }
        if (biosig[0][key] == "unknown effect"){
            biosig[1][key] = biosig[1][key] + 3;
        }
    }
    console.log(biosig);
    console.log(desired);
    renderBiosig(cmpdName,"",biosig,desired);
}
