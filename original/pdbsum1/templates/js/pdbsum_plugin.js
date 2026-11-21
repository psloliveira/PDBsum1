var metals_colours =
    {"green":"#00FF00", "cyan":"#00FFFF", "orange":"#FFA500", 
     "pink":"#FF69B4", "brown":"#8B4513"};

var colour_dict =
    {'A':'#cc66ff', 'B':'#ff0000', 'C':'#cc9900', 'D':'#ff69b4', 'E':'#0000ff', 
     'F':'#00ff00', 'G':'#33ccff', 'H':'#ffff00', 'I':'#cc66ff', 'J':'#ff0000', 
     'K':'#cc9900', 'L':'#ff69b4', 'M':'#0000ff', 'N':'#00ff00', 'O':'#33ccff', 
     'P':'#ffff00', 'Q':'#cc66ff', 'R':'#ff0000', 'S':'#cc9900', 'T':'#ff69b4', 
     'U':'#0000ff', 'V':'#00ff00', 'W':'#33ccff', 'X':'#ffff00', 'Y':'#cc66ff', 
     'Z':'#ff0000', '1':'#cc9900', '2':'#ff69b4', '3':'#0000ff', '4':'#00ff00', 
     '5':'#33ccff', '6':'#ffff00', ' ':'#cc66ff', '7':'#ff0000', '8':'#cc9900', 
     '9':'#ff69b4', '0':'#0000ff', '!':'#00ff00', '#':'#33ccff', '_':'#ffff00', 
     '-':'#cc66ff', '=':'#ff0000', ':':'#cc9900', '<':'#ff69b4', '>':'#0000ff', 
     '|':'#00ff00', ' ':'#cc66ff'};

var glviewer = null;
var apply_styles = function(viewer, dfn_file){
        
    var chains_prot = [];
    var chains_calpha = [];
    var chains_nuc = [];
    var ligands_res = [];
    var metals_res = [];

    var dfn_f = dfn_file;
    var token = dfn_f.split("\n");

    for (var i = 0; i < token.length; i++) {
        if (token[i][0] === "P"){                
            chains_prot.push(token[i][1]);
        } else if (token[i][0] === "C"){
            chains_calpha.push(token[i][1]);
        } else if (token[i][0] === "N"){
            chains_nuc.push(token[i][1]);
        } else if (token[i][0] === "L"){

            var charExp = /\D*/;
            var numbExp = /\d+\D*/g;
            all_resids = token[i].split(' ');
           
            for (var j = 1; j < all_resids.length; j++){
                if(all_resids[j].indexOf(':') < 0){
                    var chain = ' ';
                    var res_namenum = all_resids[j];
                } else {
                    var chain = all_resids[j].split(':')[1][0];
                    var res_namenum = all_resids[j].split(':')[0];
                }        
                                
                if (res_namenum[0] === "[") {
                    var tmp = res_namenum.split(']')[0];
                    var resname = tmp.split('[')[1];
                    var resnum = res_namenum.split(']')[1];
                } else {
                    var resname = res_namenum.match(charExp)[0];
                    var resnum = res_namenum.match(numbExp)[0];
                }
                ligands_res.push([resname,resnum,chain]);
            }
                        
        } else if (token[i][0] === "M"){
            if (token[i][1] === " "){
                metals_res.push([token[i].slice(1, token.length).split(" ")[1],
                                 token[i].split(" ")[2]]);
            } else {
                metals_res.push([token[i].slice(1, token.length).split(" ")[0],
                                 token[i].split(" ")[1]]);
            }
        }
    }

    for (var i = 0; i < chains_prot.length; i++) {

        // Colour by chain
        glviewer.setStyle({chain:chains_prot[i]}, 
                          {cartoon:{color:colour_dict[chains_prot[i]]}});
     }

    for (var i = 0; i < chains_nuc.length; i++) {
        glviewer.setStyle({chain:chains_nuc[i]}, 
                          {cartoon:{color:colour_dict[chains_nuc[i]]}}); // nucleic acid
    }
        
    for (var i = 0; i < ligands_res.length; i++) {
        glviewer.setStyle({chain:ligands_res[i][2], resi:ligands_res[i][1], resn:ligands_res[i][0]}, 
			  {sphere:{}}); // ligands
        glviewer.setStyle({chain:ligands_res[i][2], resi:ligands_res[i][1], resn:ligands_res[i][0], elem:'C'}, 
			  {sphere:{color:'#333333'}});
        glviewer.setStyle({chain:ligands_res[i][2], resi:ligands_res[i][1], resn:ligands_res[i][0], elem:'H'}, 
			  {sphere:{hidden:true}});
    }

    for (var i = 0; i < metals_res.length; i++) {
        glviewer.setStyle({resn:metals_res[i][0]}, 
                          {sphere:{color:metals_colours[metals_res[i][1]]}});
    }
}

$(document).ready(function() {
    
    var scripts = document.getElementsByTagName('script');

    for (var i = 0; i < scripts.length; i++) {
        if (scripts[i].accessKey === 'pdbsum_js'){
            var script = scripts[i];
        }
    } 
    
    var pdb_id = script.getAttribute('pdb_code');    

    // Read in the rasmol.dfn definitions
    dfn_file =
	//#insert dfn file here
    "\n";

    // Read in the coordinates from the PDB file
    pdb_coords =
	//#insert PDB file here
    "\n";

    var element = $('#gldiv');
    glviewer = $3Dmol.createViewer(element);
    
    m = glviewer.addModel(pdb_coords, "pdb");
    
    var m = glviewer.getModel();
    var atoms = m.selectedAtoms({});

    // Set background to white
    glviewer.setBackgroundColor(0xffffff); 

    // Display the structure
    apply_styles(glviewer, dfn_file);
    glviewer.rotate(180,"x");
    glviewer.zoomTo();
    glviewer.render();
    glviewer.zoom(1.2, 1000);
});

