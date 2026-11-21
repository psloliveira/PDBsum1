// Define metal colours
var metals_colours =
    {"green":"#00FF00", "cyan":"#00FFFF", "orange":"#FFA500", 
     "pink":"#FF69B4", "brown":"#8B4513"}

// Initialise
var glviewer = null;
var natoms = 0;

//
// Function to count number of atoms in the PDB file
//
var count_atoms = function(pdb_coords)
{
    var natoms = 0;

    var pdb_f = pdb_coords;
    var token = pdb_f.split("\n");

    //
    // Loop over the lines in the file
    //
    for (var i = 0; i < token.length; i++)
    {
        //
        // Get the record type
        //
        var rec_type = token[i].substring(0,6).trim();

        //
        // If first field is ATOM or HETATM, increment count
        //
        if (rec_type === "ATOM" || rec_type === "HETATM")
        {
            natoms++;
        }
    }

    return natoms
}

//
// Function to count number of residue in the PDB file
//
var count_residues = function(pdb_coords)
{
    var nres = 0;

    var pdb_f = pdb_coords;
    var token = pdb_f.split("\n");

    var last_chain = "";
    var last_res_name = "";
    var last_res_num = "";

    //
    // Loop over the lines in the file
    //
    for (var i = 0; i < token.length; i++)
    {
        //
        // Extract residue details
        //
        var rec_type = token[i].substring(0,6).trim();
        var res_name = token[i].substring(17,20).trim();
        var chain = token[i].substring(21,22);
        var res_num = token[i].substring(22,27).trim();

        // If residue same, ignore
        if (chain === last_chain && 
            res_name === last_res_name &&
            res_num === last_res_num)
        {
        }

        // Otherwise, increment residue count
        else if (rec_type === 'ATOM' || rec_type === 'HETATM')
        {
            nres++;
// debug
//            console.log(nres + '. ' + chain + res_name + res_num);
// debug
        }

        // Save details
        last_chain = chain;
        last_res_name = res_name;
        last_res_num = res_num;
    }

    return nres
}

//
// Add atom and residue labels
//
var add_labels = function(glviewer,pdb_coords,nres)
{
    var pdb_f = pdb_coords;
    var token = pdb_f.split("\n");

    var last_chain = "";
    var last_res_name = "";
    var last_res_num = "";

    //
    // Loop over the lines in the file
    //
    for (var i = 0; i < token.length; i++)
    {
        //
        // Initialise font size and colour
        //
        var fsize = 16;
        var fcolour = "black";

        //
        // Extract atom name, residue name, residue number and chain id
        //
        var rec_type = token[i].substring(0,6).trim();
        var atom_name = token[i].substring(12,16).trim();
        var res_name = token[i].substring(17,20).trim();
        var chain = token[i].substring(21,22);
        var res_num = token[i].substring(22,27).trim();

        // If residue same, use atom name for label
        if (nres < 2 || (chain === last_chain && 
            res_name === last_res_name &&
            res_num === last_res_num))
        {
            var labelText = '  ' + atom_name;
        }

        // Otherwise add residue to atom name
        else
        {
            var labelText = '  ' + atom_name + ' - ' + res_name
            + res_num + '(' + chain + ')';;
            fsize = 20;
            fcolour = "blue";
        }

        // Define atom to which label is to be attached
        var selector = {chain: chain, resi: res_num, atom: atom_name};
        
        // Define label appearance
        var labelSpec = {fontSize: fsize, showBackground: false,
                         fontColor: fcolour};

        // Display label if have a valid record
        if (rec_type === 'ATOM' || rec_type === 'HETATM')
        {
          var label_id = glviewer.addLabel(labelText, labelSpec, selector);
        }

        // Save details
        last_chain = chain;
        last_res_name = res_name;
        last_res_num = res_num;
    }
}

//
// Display document
//
$(document).ready
(
    function()
    {
        //Initialise
	natoms = 0;

        //
        // Get the parameters with which this script was called
        //
        var scripts = document.getElementsByTagName('script');
        for (var i = 0; i < scripts.length; i++)
        {
            if (scripts[i].accessKey === 'ligand_js')
            {
                var script = scripts[i];
            }
        } 
    
        //
        // Extract the PDB code and other information
        //
        var pdb_code = script.getAttribute('pdb_code');    
        var ligand_seq = script.getAttribute('ligand_seq');
        var clig_no = script.getAttribute('clig_no');
        var labels_on = script.getAttribute('labels');

	//
        // Get the PDB coordinates
        //
        pdb_coords =
      	//#insert PDB file here
    "\n";

        //
        // Check that file not empty
        //
        natoms = count_atoms(pdb_coords);
// debug
            console.log('PDB contents: ' + pdb_coords);
            console.log('natoms: ' + natoms);
// debug

        //
        // Display the structure if there is anything to display
        //
        if (natoms > 0)
        {
            //
            // Count number of residues
            //
            nres = count_residues(pdb_coords);
// debug
            console.log('nres: ' + nres);
// debug

            //
            // Define the viewing window
            //
            var element = $('#gldiv');
            glviewer = $3Dmol.createViewer(element);

            //
            // Add PDB coords
            //
            m = glviewer.addModel(pdb_coords, "pdb");
    
            var m = glviewer.getModel();
            var atoms = m.selectedAtoms({});

            //
            // Set background to white
            //
            glviewer.setBackgroundColor(0xffffff); 

            //
            // If required, add atom and residue labels
            //
            if (labels_on === "TRUE")
            {
              add_labels(glviewer,pdb_coords,nres);
            }

            //
            // Show as sticks
            //
            glviewer.setStyle({model:0}, {stick:{}});
            glviewer.rotate(180,"x");
            glviewer.zoomTo();
            glviewer.render();
            glviewer.zoom(1.1, 1000);
        }
    }
);

