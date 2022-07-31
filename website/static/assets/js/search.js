var resultsList = document.getElementById('results');
var response = null;
var query;
$( "#button" ).click(function() {
resultsList.innerHTML = "";
query = document.getElementById('search').value;
response = null;
load(query);
console.log("searching...");
});

function capitalizeFirstLetter(string) {
  return string.charAt(0).toUpperCase() + string.slice(1);
}

function search(){
    var term = document.getElementById('search').value;
    queryJSON = {
        "term": term
    };
    $.ajax({
        url: '/search',
        type: 'POST',
        contentType: 'application/json; charset=utf-8',
        datatype: "json",
        data: JSON.stringify(queryJSON),
        success: function (data) {
            console.log("Recieving search results: " + JSON.stringify(data));
            response = data;
            renderResults(response);
            //document.getElementById("loading").remove();
        },
        error: function () {
            console.log('Error');
        }
    });
}

const timer = ms => new Promise(res => setTimeout(res, ms))
async function load(query){
    search();
    var loadingContent = document.createElement('div');
    loadingContent.innerHTML = '<p class="loading" id="loading">Loading.</p>'
    resultsList.prepend(loadingContent);
    var inc = 0;
    while(response == null){
        setTimeout(function(){
        if (response == null){
            loadingContent.innerHTML = '<p class="loading" id="loading">Loading..</p>';
        } else {
            document.getElementById("loading").remove();
        }
        }, 1000);
        setTimeout(function(){
        if (response == null){
            loadingContent.innerHTML = '<p class="loading" id="loading">Loading...</p>';
        } else {
            document.getElementById("loading").remove();
        }
        }, 2000);
        setTimeout(function(){
        if (response == null){
            loadingContent.innerHTML = '<p class="loading" id="loading">Loading</p>';
        } else {
            document.getElementById("loading").remove();
        }
        }, 3000);
        setTimeout(function(){
        if (response == null){
            loadingContent.innerHTML = '<p class="loading" id="loading">Loading.</p>';
        } else {
            document.getElementById("loading").remove();
        }
        }, 4000);
        await timer(4000);
        if (response != null){
            document.getElementById("loading").remove();
        }
    }
}


function displayToggle(cmpdName,toggle){
    var props = document.getElementById(cmpdName+"-props");
    var biosig = document.getElementById(cmpdName+"-biosig");
    var assay = document.getElementById(cmpdName+"-assays");
    var synth = document.getElementById(cmpdName+"-synth");
    if(toggle == "props"){
        props.style.display = "block";
        biosig.style.display = "none";
        assay.style.display = "none";
        synth.style.display = "none";
    }
    if(toggle == "biosig"){
        props.style.display = "none";
        biosig.style.display = "block";
        assay.style.display = "none";
        synth.style.display = "none";
    }
    if(toggle == "assay"){
        props.style.display = "none";
        biosig.style.display = "none";
        assay.style.display = "block";
        synth.style.display = "none";
    }
    if(toggle == "synth"){
        props.style.display = "none";
        biosig.style.display = "none";
        assay.style.display = "none";
        synth.style.display = "block";
    }
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
              name: 'Partial agonism (0-1)',
              color: '#00A100'
            },
            {
              from: 1,
              to: 2,
              name: 'Agonism (1-2)',
              color: '#128FD9'
            },
            {
              from: 2,
              to: 3,
              name: 'Antagonism (2-3)',
              color: '#FFB200'
            },
            {
              from: 3,
              to: 4,
              name: 'Unknown effect (3-4)',
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
  tickAmount: 0,
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

    var chart = new ApexCharts(document.querySelector("#"+cmpdName+"-biosig"), options);
    chart.render();
}


function renderResults(results){
    var numEntries = results.ids.length;
    let entries = results.ids;
    let props = results.props;
    let assays = results.assays;
    let biosig = results.biosigs;
    let synths = results.synths;
    let desired = results.desired;
    var res = document.getElementById('results');
    for(let i = 0; i < numEntries; i++){
        // LATER DO NOT INDEX INTO PROPS LIST; KEY IN BY INCHI KEY
        let cmpdName = capitalizeFirstLetter(entries[i].name);
        let inchikey = entries[i]["_id"];
        let smiles = entries[i].smiles;
        let MolecularFormula = props[i].pubchem.MolecularFormula;
        let inchi = entries[i].inchi;
        let molwt = props[i].pubchem.ExactMass;
        let hbonddonors = props[i].pubchem.HBondDonorCount;
        let hbondacceptors = props[i].pubchem.HBondAcceptorCount;
        let logp = props[i].rdkit.MolLogP;
        let assaysJSON = JSON.stringify(assays,null,"\t");
        let image = props[i].svg;
        let assaysList = "";
        for(let j=0;j<assays[i].length;j++){
            assaysList += '<p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Assay  '+assays[i][j]._id+'<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 0px;">Receptor: '+assays[i][j].receptor+'<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">Source: '+assays[i][j].source+'<br></p>';
        }
        if (assays[i].length == 0) {
            assaysList = '<p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">No assay information available<br></p>';
        }
        let synthsList = "";
        //for (let j=0; j<synths[i].length; j++){
        for (let j=0; j<Math.min(5,synths[i].length); j++){
            synthsList += '<img src="data:image/png;base64,'+synths[i][j]+'">'
        }
        if (synths[i].length == 0) {
            synthsList = '<p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">No synthetic information available<br></p>';
        }
        let entry = document.createElement("div");
        entry.innerHTML = '<div class="row" id="'+cmpdName+'" style="border-radius: 25px;border-color: var(--bs-indigo);"><div class="col-3 d-xl-flex justify-content-xl-end"><!--img width="200px" height="200px" src="data:image/svg;base64,"-->'+image+'</div><div class="col-9"><a style="color: grey;font-size: 30px;text-align: left;margin-bottom: 5px;" href="/summary/'+inchikey+'">'+cmpdName+'</a><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">SMILES<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+smiles+'<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">InChI Key<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+inchikey+'<br></p><div class="container"><div class="row"><div class="col-md-3"><button class="btn btn-primary border-white collapsed" type="button" style="width: 100%;background: var(--bs-gray-800);border-radius: 0px;border-style: none;" data-bs-toggle="collapse" data-bs-target="#'+cmpdName+'-collapse" aria-expanded="false" aria-controls="'+cmpdName+'-collapse" onclick="displayToggle(&#39;'+cmpdName+'&#39;,&#39;props&#39;)">Properties</button></div><div class="col-md-3"><button class="btn btn-primary border-white collapsed" type="button" style="width: 100%;background: var(--bs-gray-800);border-radius: 0px;border-style: none;" data-bs-toggle="collapse" data-bs-target="#'+cmpdName+'-collapse" aria-expanded="false" aria-controls="'+cmpdName+'-collapse" onclick="displayToggle(&#39;'+cmpdName+'&#39;,&#39;biosig&#39;)">Biosignature</button></div><div class="col-md-3"><button class="btn btn-primary border-white collapsed" type="button" style="width: 100%;background: var(--bs-gray-800);border-radius: 0px;border-style: none;" data-bs-toggle="collapse" data-bs-target="#'+cmpdName+'-collapse" aria-expanded="false" aria-controls="'+cmpdName+'-collapse" onclick="displayToggle(&#39;'+cmpdName+'&#39;,&#39;assay&#39;)">Assays</button></div><div class="col-md-3"><button class="btn btn-primary border-white collapsed" type="button" style="width: 100%;background: var(--bs-gray-800);border-radius: 0px;border-style: none;" data-bs-toggle="collapse" data-bs-target="#'+cmpdName+'-collapse" aria-expanded="false" aria-controls="'+cmpdName+'-collapse" onclick="displayToggle(&#39;'+cmpdName+'&#39;,&#39;synth&#39;)">Synthesis</button></div></div></div><div class="row" style="border-style: none;"><div class="col" style="border-style: none;"><div class="row"><div class="col"><div id="'+cmpdName+'-collapse" class="collapse" style="padding: 10px;"><div id="'+cmpdName+'-biosig" style="min-height: 365px; display: none;"><div id="apexchartsqfpbc9xf" class="apexcharts-canvas apexchartsqfpbc9xf" style="width: 0px; height: 350px;"><svg id="SvgjsSvg1181" width="0" height="350" xmlns="http://www.w3.org/2000/svg" version="1.1" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:svgjs="http://svgjs.dev" class="apexcharts-svg" xmlns:data="ApexChartsNS" transform="translate(0, 0)" style="background: transparent;"><g id="SvgjsG1184" class="apexcharts-annotations"></g><g id="SvgjsG1183" class="apexcharts-inner apexcharts-graphical"><defs id="SvgjsDefs1182"></defs></g></svg><div class="apexcharts-legend"></div></div></div><div id="'+cmpdName+'-props" style="display: block;"><div class="row"><div class="col"><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">InChI<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+inchi+'<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Molecular Formula<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+MolecularFormula+'<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">CAS<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+'Not returned'+'<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Molecular Weight<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+molwt+'<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Hydrogen Bond Donor Count<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+hbonddonors+'<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Hydrogen Bond Acceptor Count<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+hbondacceptors+'<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">logP<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+logp+'<br></p></div></div></div><div id="'+cmpdName+'-assays"><div class="row"><div class="col">'+assaysList+'<!--pre><code class="prettyprint">'+assaysJSON+'</code></pre--></div></div></div><div id="'+cmpdName+'-synth">'+synthsList+'</div></div></div></div></div></div></div><div class="col"><hr></div></div>'
        res.appendChild(entry);
        
        for (var key of Object.keys(desired[0])){
            if (desired[1][key] == -1){
                continue;
            }
            if (desired[0][key] == "partial agonist"){
                desired[0][key] = desired[0][key];
            }
            if (desired[0][key] == "binder"){
                desired[1][key] = desired[1][key]  + 3;
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
        for (var key of Object.keys(biosig[i][0])){
            if (desired[1][key] == -1){
                continue;
            }
            if (biosig[i][key] == "partial agonist"){
                biosig[i][key] = biosig[i][key];
            }
            if (biosig[i][0][key] == "binder"){
                biosig[i][1][key] = biosig[i][1][key] + 3;
            }
            if (biosig[i][0][key] == "agonism"){
                biosig[i][1][key] = biosig[i][1][key] + 1;
            }
            if (biosig[i][0][key] == "antagonism"){
                biosig[i][1][key] = biosig[i][1][key] + 2;
            }
            if (biosig[i][0][key] == "unknown effect"){
                biosig[i][1][key] = biosig[i][1][key] + 3;
            }
        }
        console.log(biosig[i]);
        console.log(desired);
        renderBiosig(cmpdName,"",biosig[i],desired);

    }
    
}


