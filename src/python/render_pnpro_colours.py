# pnpro_renderer.py
import xml.etree.ElementTree as ET
import copy
import logging

def parse_xml_string(xml_string):
    return ET.fromstring(xml_string)

def write_xml(root, filename="output.PNPRO"):
    tree = ET.ElementTree(root)
    def indent(elem, level=0):
        i = "\n" + level * "  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for child in elem:
                indent(child, level + 1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i
    indent(root)
    tree.write(filename, encoding="UTF-8", xml_declaration=True)
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()
    content = content.replace('<?xml version="1.0" encoding="UTF-8"?>',
                             '<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
    if '<!-- This project file has been saved by the New GreatSPN Editor' not in content:
        content = content.replace('<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
                                 '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n<!-- This project file has been saved by the New GreatSPN Editor, v.100 -->')
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(content)

def identify_organisms(nodes):
    organisms = set()
    for node in nodes:
        name = node.get('name', '')
        if name.startswith('n_') or name.startswith('biomass_e_'):
            parts = name.split('_')
            if len(parts) > 1:
                organisms.add(parts[-1])
        elif '_in_' in name or '_out_' in name:
            parts = name.split('_')
            if len(parts) > 1:
                organisms.add(parts[-1])
        elif name.startswith(('Dup_', 'Death_', 'Starv_')):
            parts = name.split('_')
            if len(parts) > 1:
                organisms.add(parts[-1])
    return sorted(list(organisms))

def identify_boundary_metabolites(nodes):
    metabolites = set()
    for node in nodes:
        if node.tag == 'place':
            name = node.get('name', '')
            if not name.startswith('n_') and not name.startswith('biomass_e_'):
                metabolites.add(name)
    return sorted(list(metabolites))

def create_template_with_organisms(blank_root, organisms):
    # We're going to use the blank template exactly as is
    # This avoids issues with modifying template nodes/arcs
    return copy.deepcopy(blank_root)

def add_arc(edges, head, tail, kind="INPUT", mult=None):
    """Helper to add a new arc with proper attributes"""
    arc = ET.SubElement(edges, 'arc')
    arc.set('head', head)
    arc.set('tail', tail)
    arc.set('kind', kind)
    if mult:
        arc.set('mult', str(mult))
    return arc

def create_petri_net_from_scratch(blank_root, input_root, organisms, boundary_mets):
    """Creates a new Petri net using the input nodes but with a clean structure"""
    template = copy.deepcopy(blank_root)
    
    # Clear all existing nodes and edges to start fresh
    gspn = template.find('.//gspn')
    nodes_elem = gspn.find('./nodes')
    if nodes_elem is not None:
        for child in list(nodes_elem):
            nodes_elem.remove(child)
    
    edges_elem = gspn.find('./edges')
    if edges_elem is not None:
        for child in list(edges_elem):
            edges_elem.remove(child)
    else:
        edges_elem = ET.SubElement(gspn, 'edges')
    
    # Create text boxes for organization
    for i, org in enumerate(organisms):
        side = 'left' if i % 2 == 0 else 'right'
        
        # Add organizational boxes
        community_box = ET.SubElement(nodes_elem, 'text-box')
        community_box.set('name', f'community_box_{side}')
        community_box.set('border-color', '#000000')
        community_box.set('fill-color', '#f3f6ff')
        community_box.set('text-color', '#000000')
        community_box.set('width', '30.0')
        community_box.set('height', '20.0')
        community_box.set('x', '2.0' if side == 'left' else '72.0')
        community_box.set('y', '1.3888890000000007')
        
        # FBA wrapper box
        fba_box = ET.SubElement(nodes_elem, 'text-box')
        fba_box.set('name', f'{org}_{side}_fba_wrappers_box')
        fba_box.set('border-color', '#000000')
        fba_box.set('fill-color', '#fffcf3')
        fba_box.set('text-color', '#000000')
        fba_box.set('width', '25.0')
        fba_box.set('height', '15.0')
        fba_box.set('x', '4.5' if side == 'left' else '74.5')
        fba_box.set('y', '3.8888890000000007')
        
        # Module box
        module_box = ET.SubElement(nodes_elem, 'text-box')
        module_box.set('name', f'{org}_{side}_module_box')
        module_box.set('border-color', '#000000')
        module_box.set('fill-color', '#f9f3ff' if side == 'left' else '#f4fff3')
        module_box.set('text-color', '#000000')
        module_box.set('width', '6.0')
        module_box.set('height', '12.0')
        module_box.set('x', '10.0' if side == 'left' else '80.0')
        module_box.set('y', '5.388888999999999')
        
        # Boundary transitions box
        boundary_box = ET.SubElement(nodes_elem, 'text-box')
        boundary_box.set('name', f'{org}_{side}_boundary_transitions_box')
        boundary_box.set('border-color', '#000000')
        boundary_box.set('fill-color', '#f4fff3' if side == 'left' else '#f9f3ff')
        boundary_box.set('text-color', '#000000')
        boundary_box.set('width', '6.0')
        boundary_box.set('height', '12.0')
        boundary_box.set('x', '18.0' if side == 'left' else '88.0')
        boundary_box.set('y', '5.388888999999999')
    
    # Add boundary metabolites box
    boundary_met_box = ET.SubElement(nodes_elem, 'text-box')
    boundary_met_box.set('name', 'boundary_metabolites_box')
    boundary_met_box.set('border-color', '#000000')
    boundary_met_box.set('fill-color', '#fff3fc')
    boundary_met_box.set('text-color', '#000000')
    boundary_met_box.set('width', '5.0')
    boundary_met_box.set('height', str(max(10.0, len(boundary_mets) * 2.5)))
    boundary_met_box.set('x', '49.5')
    boundary_met_box.set('y', '6.388888999999999')
    
    # Add all places, transitions and edges from input file
    # but preserve our clean structure and avoid incomplete references
    input_nodes = input_root.find('.//gspn/nodes')
    input_edges = input_root.find('.//gspn/edges')
    
    # First add all nodes with proper positioning
    node_map = {}  # To track original -> new node mapping
    
    # Define base positions for different node types
    base_positions = {
        'n_': {'x': 5.0, 'y': 7.0, 'dx': 8.0, 'dy': 2.0},
        'biomass_e_': {'x': 5.0, 'y': 9.0, 'dx': 8.0, 'dy': 2.0},
        'glc__D_e': {'x': 49.5, 'y': 7.0, 'dx': 0, 'dy': 2.0},
        'lcts_e': {'x': 49.5, 'y': 9.0, 'dx': 0, 'dy': 2.0},
        'Dup_': {'x': 6.0, 'y': 12.0, 'dx': 8.0, 'dy': 2.0},
        'Death_': {'x': 10.0, 'y': 12.0, 'dx': 8.0, 'dy': 2.0},
        'Starv_': {'x': 8.0, 'y': 14.0, 'dx': 8.0, 'dy': 2.0},
        'EX_biomass_e_in_': {'x': 16.0, 'y': 8.0, 'dx': 8.0, 'dy': 2.0},
        'EX_biomass_e_out_': {'x': 20.0, 'y': 8.0, 'dx': 8.0, 'dy': 2.0},
        'EX_glc__D_e_in_': {'x': 16.0, 'y': 12.0, 'dx': 8.0, 'dy': 2.0},
        'EX_glc__D_e_out_': {'x': 20.0, 'y': 12.0, 'dx': 8.0, 'dy': 2.0},
        'EX_lcts_e_in_': {'x': 16.0, 'y': 16.0, 'dx': 8.0, 'dy': 2.0},
        'EX_lcts_e_out_': {'x': 20.0, 'y': 16.0, 'dx': 8.0, 'dy': 2.0},
    }
    
    # Track indices for each node type
    node_indices = {prefix: 0 for prefix in base_positions.keys()}
    
    # Add places first
    for place in input_nodes.findall('./place'):
        new_place = copy.deepcopy(place)
        name = place.get('name', '')
        
        # Find appropriate position based on node type
        for prefix, position in base_positions.items():
            if name.startswith(prefix):
                idx = node_indices[prefix]
                new_place.set('x', str(position['x'] + idx * position['dx']))
                new_place.set('y', str(position['y'] + idx * position['dy']))
                node_indices[prefix] += 1
                break
        
        nodes_elem.append(new_place)
        node_map[name] = new_place
    
    # Add transitions
    for transition in input_nodes.findall('./transition'):
        new_transition = copy.deepcopy(transition)
        name = transition.get('name', '')
        
        # Find appropriate position based on node type
        for prefix, position in base_positions.items():
            if name.startswith(prefix):
                idx = node_indices[prefix]
                new_transition.set('x', str(position['x'] + idx * position['dx']))
                new_transition.set('y', str(position['y'] + idx * position['dy']))
                node_indices[prefix] += 1
                break
        
        nodes_elem.append(new_transition)
        node_map[name] = new_transition
    
    # Now add edges, but only if both source and target exist
    if input_edges is not None:
        for arc in input_edges:
            head = arc.get('head', '')
            tail = arc.get('tail', '')
            kind = arc.get('kind', 'INPUT')
            mult = arc.get('mult')
            
            # Only create arcs between nodes that exist in our new structure
            if head in node_map and tail in node_map:
                add_arc(edges_elem, head, tail, kind, mult)
    
    return template

def render_pnpro_layout(blank_content, input_content, output_filename="output.PNPRO"):
    try:
        blank_root = parse_xml_string(blank_content)
        input_root = parse_xml_string(input_content)
        input_nodes = input_root.findall('.//gspn/nodes/*')

        logging.info("Found %d nodes in input file.", len(input_nodes))

        organisms = identify_organisms(input_nodes)
        boundary_mets = identify_boundary_metabolites(input_nodes)

        logging.info("Organisms: %s", organisms)
        logging.info("Boundary metabolites: %s", boundary_mets)

        # Complete redesign of approach - create a clean Petri net from scratch
        # instead of trying to adapt the template
        final_pn = create_petri_net_from_scratch(blank_root, input_root, organisms, boundary_mets)
        
        write_xml(final_pn, output_filename)
        return f"Successfully generated PNPRO layout: {output_filename}"

    except Exception as e:
        logging.exception("Failed to generate PNPRO layout")
        return f"Error generating PNPRO layout: {str(e)}"

if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO)
    
    if len(sys.argv) < 3:
        print("Usage: python pnpro_renderer.py <blank_file> <input_file> [output_file]")
        sys.exit(1)

    blank_file = sys.argv[1]
    input_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else "output.PNPRO"

    with open(blank_file, 'r', encoding='utf-8') as f:
        blank_content = f.read()

    with open(input_file, 'r', encoding='utf-8') as f:
        input_content = f.read()

    result = render_pnpro_layout(blank_content, input_content, output_file)
    print(result)