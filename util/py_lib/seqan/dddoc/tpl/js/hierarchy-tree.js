// Namespace seqan.
var seqan;
if (seqan === undefined)
    seqan = {};
// Namespace seqan.doc.
if (seqan.doc === undefined)
    seqan.doc = {};

// Constants.
seqan.doc.constants = {
    FONT: 'Arial',
    FONT_SIZE: 12,
    FONT_SIZE_ROOT: 14,
    TEXT_PADDING: 6,
    NODE_SPACING: 20,
    MARGIN: 10.5,
    LIMIT_FIRST: 5,
    LIMIT_AFTER: 2,
    ARROW_DELTA: 4,

    SCHEMES: {
        'concept': {
            FILL_STYLE: '#DDF',
            STROKE_STYLE: '2px #666 dashed',
            TEXT_COLOR: '#333',
        },
        'class': {
            FONT_STYLE: 'bolder italic',
            FILL_STYLE: '#FFC',
            STROKE_STYLE: '2px #666',
            TEXT_COLOR: '#000000'
        },
        'spec': {
            FONT_STYLE: '',
            FILL_STYLE: '#FFC',
            STROKE_STYLE: '2px #666',
            TEXT_COLOR: '#000000'
        },
    }
};

/*
  Represents one entry in the class hierarchy.

  Note that only the root node should have both parents and children.
*/
seqan.doc.HierarchyNode = function(title, level, type, link)
{
    this.level = level;
    this.title = title;
    this.parents = [];
    this.children = [];
    this.dim = null;
    this.label_dim = null;
    this.subtree_dim = null;
    // Dimension of the trail up from here, if root node. [0, 0] otherwise.
    this.uptrail_dim = [0,0];
    this.pos = null;
    this.label_pos = null;
    this.row_widths = null;
    this.row_heights = null;
    this.type = type;
    this.link = link;
    this.isMouseOver = false;
}

/*
  Set the dimensions of a HierarchyNode.

  The dimensions are computed from the label and subtree dimensions.
*/
seqan.doc.HierarchyNode.prototype.updateDim = function()
{
    var c = seqan.doc.constants;
    var w = this.label_dim[0];
    var h = this.label_dim[1];
    if (this.subtree_dim[0] > 0)
    {
        w = Math.max(w, this.subtree_dim[0]);
        h = h + this.subtree_dim[1];
    }
    if (this.uptrail_dim[0] > 0)
    {
        w = Math.max(w, this.uptrail_dim[0]);
        h = h + this.uptrail_dim[1];
    }
    this.dim = [w, h];
}

/*
  Return true if the given coordinates are within the shape.
*/
seqan.doc.HierarchyNode.prototype.isPointInLabel = function(x, y)
{
    x -= seqan.doc.constants.MARGIN;
    y -= seqan.doc.constants.MARGIN;
    return (x >= this.label_pos[0] && x <= this.label_pos[0] + this.label_dim[0] &&
            y >= this.label_pos[1] && y <= this.label_pos[1] + this.label_dim[1]);
}


/*
  Convert a JSON/object representation of HierarchyNode objects into HierarchyNode object tree.
*/
seqan.doc.transformHierarchy = function(obj, level)
{
    if (typeof(level) == 'undefined') level = 0;

    var type = (typeof(obj['type']) == 'undefined') ? 'spec' : obj.type;
    
    var res = new seqan.doc.HierarchyNode(obj['title'], level, type, obj['link']);
    if (typeof(obj['children']) != 'undefined')
        for (k in obj['children'])
            res.children.push(seqan.doc.transformHierarchy(obj['children'][k], level + 1));
    if (typeof(obj['parents']) != 'undefined')
        for (k in obj['parents'])
            res.parents.push(seqan.doc.transformHierarchy(obj['parents'][k], level + 1));
    return res;
}

// Code taken from http://stackoverflow.com/questions/1134586
seqan.doc.measureTextHeight = function(ctx, left, top, width, height)
{

    // Draw the text in the specified area
    ctx.save();
    ctx.translate(left, top + Math.round(height * 0.8));
    ctx.fillText('gM', 0, 0); // This seems like tall text...  Doesn't it?
    ctx.restore();

    // Get the pixel data from the canvas
    var data = ctx.getImageData(left, top, width, height).data,
        first = false, 
        last = false
        r = height,
        c = 0;

    // Find the last line with a non-white pixel
    while(!last && r) {
        r--;
        for(c = 0; c < width; c++) {
                if(data[r * width * 4 + c * 4 + 3]) {
                        last = r;
                        break;
                }
        }
    }

    // Find the first line with a non-white pixel
    while(r) {
        r--;
        for(c = 0; c < width; c++) {
                if(data[r * width * 4 + c * 4 + 3]) {
                        first = r;
                        break;
                }
        }

        // If we've got it then return the height
        if(first != r) return last - first;
    }

    // We screwed something up...  What do you expect from free code?
    return 0;
}

/*
  Recursively compute the label dimensions for the "upward trail" from the current node.
*/
seqan.doc.computeLayoutDimsUpTrail = function(node, ctx, isRoot)
{
    var c = seqan.doc.constants;

    // Measure text height.
    ctx.font = seqan.doc.constants.FONT;
    var fnt = seqan.doc.constants.FONT;
    if (isRoot)
        fnt = seqan.doc.constants.FONT_SIZE_ROOT + 'px ' + fnt;
    else
        fnt = seqan.doc.constants.FONT_SIZE + 'px ' + fnt;
    fnt = seqan.doc.constants.SCHEMES[node.type.toLowerCase()].FONT_STYLE + ' ' + fnt;
    ctx.font = fnt;
    var height = isRoot ? seqan.doc.constants.FONT_SIZE_ROOT : seqan.doc.constants.FONT_SIZE;
    //console.log(isRoot, ctx.font, height);

    // Compute dimensions of the label.
    var tm = ctx.measureText(node.title);
    var w = 2 * c.TEXT_PADDING + tm.width;
    var h = 2 * c.TEXT_PADDING + height;
    node.label_dim = [w, h];
    node.subtree_dim = [0, 0];

    // Recurse (linearly) up the trail.
    if (node.parents.length > 0)
    {
        seqan.doc.computeLayoutDimsUpTrail(node.parents[0], ctx, false);

        node.uptrail_dim[0] = Math.max(node.parents[0].uptrail_dim[0], node.parents[0].label_dim[0]);
        //console.log(node.uptrail_dim, node.parents[0].label_dim);
        node.uptrail_dim[1] = node.parents[0].uptrail_dim[1] + c.NODE_SPACING + node.label_dim[1];
        //console.log(node.uptrail_dim, node.parents[0].label_dim);
    }

    node.updateDim();
}

/*
  Recursively compute the label and subtree dimensions of a Hierarchy Node.

  ctx -- Canvas context to use for text measures.
 */
seqan.doc.computeLayoutDims = function(node, ctx, level/*=0*/, isRoot/*=false*/)
{
    if (typeof(level) == 'undefined') level = 0;
    if (typeof(isRoot) == 'undefined') isRoot = false;
    var c = seqan.doc.constants;

    // Measure text height.
    ctx.font = seqan.doc.constants.FONT;
    var fnt = seqan.doc.constants.FONT;
    if (isRoot)
        fnt = seqan.doc.constants.FONT_SIZE_ROOT + 'px ' + fnt;
    else
        fnt = seqan.doc.constants.FONT_SIZE + 'px ' + fnt;
    fnt = seqan.doc.constants.SCHEMES[node.type.toLowerCase()].FONT_STYLE + ' ' + fnt;
    ctx.font = fnt;
    var height = isRoot ? seqan.doc.constants.FONT_SIZE_ROOT : seqan.doc.constants.FONT_SIZE;
    //console.log(isRoot, ctx.font, height);

    // Compute dimensions of the label.
    var tm = ctx.measureText(node.title);
    var w = 2 * c.TEXT_PADDING + tm.width;
    var h = 2 * c.TEXT_PADDING + height;
    //console.log(h, c.TEXT_PADDING, height);
    node.label_dim = [w, h];

    // Base case: No more children.
    if (node.children.length == 0)
    {
        node.subtree_dim = [0, 0];
        node.updateDim();
        return 0;  // OK, no children.
    }

    // Recurse on children.
    var subtree = node.children;
    for (i in subtree)
        seqan.doc.computeLayoutDims(subtree[i], ctx, level + 1);
    // Now, compute layout in rows.  The constant LIMIT_FIRST gives the highest
    // number of elements in the first, LIMIT_AFTER the number of elements in
    // following rows.
    var ws = [0];  // Widths of the rows.
    var hs = [0];  // Heights of the rows.
    var j = 0;     // Index of the current row.
    for (i in subtree)
    {
        var child = subtree[i];
        ws[ws.length - 1] += child.dim[0];
        hs[hs.length - 1] = Math.max(hs[hs.length - 1], child.dim[1]);
        //console.log('i ==', i, 'row', Math.floor(i/c.LIMIT_FIRST), 'current height', hs[hs.length - 1])
        if ((level < 1 && j == c.LIMIT_FIRST - 1) || (level >= 1 && j == c.LIMIT_AFTER - 1))
        {
            j = 0;
            ws.push(0);
            hs.push(0);
        }
        else
        {
            j += 1;
            ws[ws.length - 1] += c.NODE_SPACING;
        }
    }
    if (hs[hs.length - 1] == 0)
    {
        hs.pop();
        ws.pop();
    }
    node.row_widths = ws;
    node.row_heights = hs;
    //console.log(ws);

    // The node's subtree dimension is now computed from the next level of trees.
    var w = Math.max.apply(Math, ws);
    var hs_sum = 0;
    for (i in hs)
        hs_sum += hs[i];
    var h = node.row_heights.length * c.NODE_SPACING + hs_sum;
    node.subtree_dim = [w, h];
    node.row_count = ws.length;
    node.updateDim();
    return 0;  // OK
}

/* Compute the layout positions, computeLayoutDims() must have been called already.
 */
seqan.doc.computeLayoutPositions = function(node, x, y, level)
{
    if (typeof(level) == 'undefined') level = 0;
    var c = seqan.doc.constants;
    node.pos = [x, y];

    // Special case: Is root, also compute positions upwards.
    function rec(node, root_width)
    {
        var delta = 0;
        if (node.parents.length > 0)
            delta = rec(node.parents[0], root_width);
        node.pos = [0, delta];
        var lx = (root_width - node.label_dim[0]) / 2;
        node.label_pos = [lx, delta];

        return delta + node.label_dim[1] + c.NODE_SPACING;
    }
    if (level == 0)
    {
        rec(node, node.dim[0]);
        //console.log(node);
        //console.log(node.title, node.label_pos)
        y = node.label_pos[1];
    }

    // Base case: Tree does not go on from here.
    if (node.children.length == 0)
    {
        node.label_pos = [x, y];
        return 0;  // OK.
    }

    // Otherwise, compute layout positions for nested trees.
    node.label_pos = [x + (node.dim[0] - node.label_dim[0]) / 2 - c.NODE_SPACING / 2, y];
    var x2 = x;
    var y2 = y + c.NODE_SPACING + node.label_dim[1];
    var x2_base = x2;
    var subtree = node.children;
    var row = 0;  // Current row.
    var j = 0;    // Index in current row.
    for (i in subtree)
    {
        var child = subtree[i];
        seqan.doc.computeLayoutPositions(child, x2, y2, level + 1);
        x2 += c.NODE_SPACING;
        x2 += child.dim[0];
        if ((level < 1 && j == c.LIMIT_FIRST - 1) || (level >= 1 && j == c.LIMIT_AFTER - 1))
        {
            //console.log("bumping y");
            j = 0;
            x2 = x2_base;
            y2 += c.NODE_SPACING + node.row_heights[row];
            row += 1
        }
        else
        {
            j += 1;
        }
    }
}

/*
  Compute the layout dimensions and positions.
*/
seqan.doc.computeLayout = function(node, ctx)
{
    seqan.doc.computeLayoutDimsUpTrail(node, ctx, true);
    seqan.doc.computeLayoutDims(node, ctx, 0, true);
    seqan.doc.computeLayoutPositions(node, 0, 0);
}

seqan.doc.drawHierarchy = function(ctx, node)
{
    var c = seqan.doc.constants;

    // Global canvas drawing configuration, recompute layout.
    ctx.font = c.FONT;
    seqan.doc.computeLayout(node, ctx);

    // Fill the background with gray.
    ctx.fillStyle = 'white';
    ctx.fillRect(0, 0, ctx.canvas.width, ctx.canvas.height);

    // Recursively render the nodes.

    // This function recursively renders one node and the line to its children.
    function renderNode(node, isRoot/*=false*/)
    {
        if (typeof(isRoot) == 'undefined') isRoot = false;
        //console.log('renderNode(', node.title, ')')
        ctx.lineCap = 'square';

        // Set font and color scheme.
        var fnt = seqan.doc.constants.FONT;
        if (isRoot)
            fnt = seqan.doc.constants.FONT_SIZE_ROOT + 'px ' + fnt;
        else
            fnt = seqan.doc.constants.FONT_SIZE + 'px ' + fnt;
        fnt = seqan.doc.constants.SCHEMES[node.type.toLowerCase()].FONT_STYLE + ' ' + fnt;
        ctx.font = fnt;
        //console.log(isRoot, ctx.font);

        // Fill box around text and fill it.
        //console.log(node);
        var colorScheme = seqan.doc.constants.SCHEMES[node.type.toLowerCase()];
        ctx.fillStyle = colorScheme.FILL_STYLE;
        ctx.strokeStyle = colorScheme.STROKE_STYLE;
        //ctx.fillRect(c.MARGIN + c.NODE_SPACING + node.label_pos[0], c.MARGIN + node.label_pos[1],
        //             node.label_dim[0], node.label_dim[1]);
        //ctx.strokeRect(c.MARGIN + c.NODE_SPACING + node.label_pos[0], c.MARGIN + node.label_pos[1],
        //               node.label_dim[0], node.label_dim[1]);
        // Draw text.
        ctx.textBaseline = 'top';
        ctx.fillStyle = colorScheme.TEXT_COLOR;
        ctx.fillText(node.title,
                     c.MARGIN + c.NODE_SPACING + c.TEXT_PADDING + node.label_pos[0],
                     c.MARGIN + c.TEXT_PADDING + node.label_pos[1]);

        // Draw bounding box.
        //ctx.strokeStyle = 'red';
        //ctx.strokeRect(c.MARGIN + c.NODE_SPACING + node.pos[0], c.MARGIN + node.pos[1], node.dim[0], node.dim[1]);

        ctx.lineWidth = 1;
        ctx.lineCap = 'round';
        if (node.children.length > 0)
        {
            // Draw line down out of text.
            var node_label_center_x = c.MARGIN + c.NODE_SPACING + Math.floor(node.label_pos[0] + node.label_dim[0] / 2.0);
            var node_label_center_y = c.MARGIN + node.label_pos[1] + node.label_dim[1];
            ctx.moveTo(node_label_center_x, node_label_center_y);
            ctx.lineTo(node_label_center_x, node_label_center_y + c.NODE_SPACING / 2.0)
            //console.log('move to', [node_label_center_x, node_label_center_y]);
            //console.log('line to', [node_label_center_x, node_label_center_y + c.NODE_SPACING / 2.0]);
            ctx.strokeStyle = '1px black';
            ctx.stroke();
            // Draw little inheritance arrow.
            ctx.moveTo(node_label_center_x, node_label_center_y);
            ctx.lineTo(node_label_center_x - c.ARROW_DELTA / 2.0, node_label_center_y + c.ARROW_DELTA);
            ctx.lineTo(node_label_center_x + c.ARROW_DELTA / 2.0, node_label_center_y + c.ARROW_DELTA);
            ctx.lineTo(node_label_center_x, node_label_center_y);
            //console.log('down');
            //console.log(node_label_center_x, node_label_center_y,
            //            node_label_center_x - c.ARROW_DELTA / 2.0, node_label_center_y + c.ARROW_DELTA,
            //            node_label_center_x + c.ARROW_DELTA / 2.0, node_label_center_y + c.ARROW_DELTA);
            ctx.fill();
            // Compute rows.
            var rows = [];
            var limit = 0;
            if (node.level< 1)
                limit = c.LIMIT_FIRST;
            else
                limit = c.LIMIT_AFTER;
            for (var i = 0; i < node.children.length; i += limit)
            {
                rows.push([i, Math.min(i + limit, node.children.length)]);
                if (i + limit >= node.children.length)
                    break;
            }
            // Draw horizontal lines and lines into children, row-wise.
            for (var i = 0; i < rows.length; i++)
            {
                var row = rows[i];
                var leftmost_x = null;
                var rightmost_x = null;
                var row_slice = node.children.slice(row[0], row[1]);
                for (j in row_slice)
                {
                    var child = row_slice[j];
                    // Compute line positions.
                    var label_center_x = c.MARGIN + c.NODE_SPACING + Math.floor(child.label_pos[0] + child.label_dim[0] / 2.0);
                    var label_center_y = c.MARGIN + child.label_pos[1];
                    // Update leftmost/rightmost x coordinates.
                    //console.log("leftmost x", leftmost_x, "rightmost x", rightmost_x);
                    if (leftmost_x === null || leftmost_x > label_center_x)
                        leftmost_x = label_center_x;
                    if (rightmost_x === null || rightmost_x < label_center_x)
                        rightmost_x = label_center_x;
                    //console.log("  leftmost x", leftmost_x, "rightmost x", rightmost_x);
                    // Draw lines.
                    ctx.moveTo(label_center_x, label_center_y);
                    ctx.lineTo(label_center_x, label_center_y - c.NODE_SPACING / 2);
                    //console.log("line", [label_center_x, label_center_y, label_center_x, label_center_y - c.NODE_SPACING / 2]);
                    ctx.stroke();

                    // Draw horizontal line.
                    if (!(leftmost_x === null || rightmost_x === null))
                    {
                        // Make sure that the line connects to the one coming from the parent.
                        if (row[0] == 0)
                        {
                            leftmost_x = Math.min(leftmost_x, node_label_center_x);
                            rightmost_x = Math.max(rightmost_x, node_label_center_x);
                        }
                        // If there is more than one row then connect to the left.
                        if (rows.length > 1)
                            leftmost_x = c.MARGIN + c.NODE_SPACING + node.pos[0] - c.NODE_SPACING + 2 * node.level;
                        // Actually draw line.
                        ctx.moveTo(leftmost_x, label_center_y - c.NODE_SPACING / 2);
                        ctx.lineTo(rightmost_x, label_center_y - c.NODE_SPACING / 2);
                        ctx.stroke()
                        //console.log("horizontal line", [leftmost_x, label_center_y - c.NODE_SPACING / 2, leftmost_x, label_center_y - c.NODE_SPACING / 2]);
                    }
                }

                // Finally, draw vertical line if more than one row.
                if (rows.length > 1)
                {
                    ctx.moveTo(c.MARGIN + c.NODE_SPACING + node.pos[0] - c.NODE_SPACING + 2 * node.level,
                               c.MARGIN + node.children[0].label_pos[1] - c.NODE_SPACING / 2);
                    ctx.lineTo(c.MARGIN + c.NODE_SPACING + node.pos[0] - c.NODE_SPACING + 2 * node.level,
                               c.MARGIN + node.children[node.children.length - 1].label_pos[1] - c.NODE_SPACING / 2);
                    ctx.stroke()
                }

                // And... stroke!
                ctx.stroke()
            }
        }
    }

    // This function renders one node and all its children.
    function recurseDown(node, isRoot/*=false*/)
    {
        if (typeof(isRoot) == 'undefined') isRoot = false;

        renderNode(node, isRoot);
        for (i in node.children)
            recurseDown(node.children[i]);
    }

    // This function renders one uptrail node.
    function renderNodeUpTrail(node, isRoot)
    {
        // Set font and color scheme.
        var fnt = seqan.doc.constants.FONT;
        fnt = seqan.doc.constants.FONT_SIZE + 'px ' + fnt;
        fnt = seqan.doc.constants.SCHEMES[node.type.toLowerCase()].FONT_STYLE + ' ' + fnt;
        ctx.font = fnt;
        //console.log(isRoot, ctx.font);

        // Fill box around text and fill it.
        //console.log(node);
        var colorScheme = seqan.doc.constants.SCHEMES[node.type.toLowerCase()];
        ctx.fillStyle = colorScheme.FILL_STYLE;
        ctx.strokeStyle = colorScheme.STROKE_STYLE;
        //ctx.fillRect(c.MARGIN + c.NODE_SPACING + node.label_pos[0], c.MARGIN + node.label_pos[1],
        //             node.label_dim[0], node.label_dim[1]);
        //ctx.strokeRect(c.MARGIN + c.NODE_SPACING + node.label_pos[0], c.MARGIN + node.label_pos[1],
        //               node.label_dim[0], node.label_dim[1]);
        // Draw text.
        ctx.textBaseline = 'top';
        ctx.fillStyle = colorScheme.TEXT_COLOR;
        ctx.fillText(node.title,
                     c.MARGIN + c.NODE_SPACING + c.TEXT_PADDING + node.label_pos[0],
                     c.MARGIN + c.TEXT_PADDING + node.label_pos[1]);
        //console.log(colorScheme)

        // Draw arrowed line to parent if any.
        ctx.lineWidth = 1;
        ctx.lineCap = 'round';
        //console.log(node.parents);
        if (node.parents.length > 0)
        {
            var node_label_center_x = c.MARGIN + c.NODE_SPACING + Math.floor(node.label_pos[0] + node.label_dim[0] / 2.0);
            var node_label_center_y = c.MARGIN + node.label_pos[1];
            ctx.moveTo(node_label_center_x, node_label_center_y);
            ctx.lineTo(node_label_center_x, node_label_center_y - c.NODE_SPACING);
            //console.log(node_label_center_x, node_label_center_y, node_label_center_x, node_label_center_y - c.NODE_SPACING);
            ctx.stroke();
            ctx.moveTo(node_label_center_x, node_label_center_y - c.NODE_SPACING);
            ctx.lineTo(node_label_center_x - c.ARROW_DELTA / 2.0, node_label_center_y - c.NODE_SPACING + c.ARROW_DELTA);
            ctx.lineTo(node_label_center_x + c.ARROW_DELTA / 2.0, node_label_center_y - c.NODE_SPACING + c.ARROW_DELTA);
            ctx.lineTo(node_label_center_x, node_label_center_y - c.NODE_SPACING);
            //console.log('line')
            //console.log(node_label_center_x, node_label_center_y - c.NODE_SPACING,
            //            node_label_center_x - c.ARROW_DELTA / 2.0, node_label_center_y - c.NODE_SPACING + c.ARROW_DELTA,
            //            node_label_center_x + c.ARROW_DELTA / 2.0, node_label_center_y - c.NODE_SPACING + c.ARROW_DELTA);
            ctx.fill();
        }
        ctx.fillStyle = 'black';
        ctx.strokeStyle = 'black';
    }

    // This function renders one node and all its parents.
    function recurseUp(node, isRoot/*=false*/)
    {
        if (typeof(isRoot) == 'undefined') isRoot = false;

        renderNodeUpTrail(node, isRoot);
        if (node.parents.length > 0)
            recurseUp(node.parents[0]);
    }
    
    recurseDown(node, true);
    recurseUp(node, true);
    return 0;  // OK
}

seqan.doc.createDivs = function(id, node)
{
    var c = seqan.doc.constants;

    function renderNode(node, isRoot/*=false*/)
    {
        //console.log(node)
        if (node.label_pos === null)
            return;  // not to be rendered (e.g. non-first concept)
        var rootClass = node.type;
        if (isRoot)
            rootClass = rootClass + ' root';
        rootClass = ' class="' + rootClass + '"';
        $(id).append('<div class="hierarchy_entry" style="position: absolute; left: ' + ($('#canvas').position().left + c.MARGIN + c.NODE_SPACING + node.label_pos[0]) + 'px; top: ' + ($('#canvas').position().top + c.MARGIN + node.label_pos[1]) + 'px; height: ' + (node.label_dim[1] - 2) + 'px; width:' + (node.label_dim[0] - 2) + '"><div' + rootClass + '><a href="' + node.link + '">' + node.title.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;") + '</a></div></div>');
    }

    function recurseDown(node, isRoot/*=false*/)
    {
        if (typeof(isRoot) == 'undefined') isRoot = false;

        renderNode(node, isRoot);
        for (i in node.children)
            recurseDown(node.children[i]);
    }

    // This function renders one node and all its parents.
    function recurseUp(node, isRoot/*=false*/)
    {
        if (typeof(isRoot) == 'undefined') isRoot = false;

        renderNode(node, isRoot);
        for (i in node.parents)
            recurseUp(node.parents[i]);
    }

    recurseDown(node, true);
    recurseUp(node, true);
}

seqan.doc.updateDivs = function(id, node)
{
    $(id).children().remove();
    seqan.doc.createDivs(id, node);
}

function initHierarchyCanvas()
{
    var canvas = document.getElementById("canvas");
    if (!canvas)
        return;  // Abort if no canvas element found.
    var ctx = canvas.getContext("2d");
    
    // Build tree in seqan.doc.HierarchyNode from data and compute the layout
    // (dimensions and positions).
    var t = seqan.doc.transformHierarchy(data);
    ctx.font = seqan.doc.constants.FONT;
    seqan.doc.computeLayout(t, ctx);
    //console.log(t);
    
    // Update canvas dimensions.
    canvas.width = t.dim[0] + 2 * seqan.doc.constants.MARGIN + seqan.doc.constants.NODE_SPACING;
    canvas.height = t.dim[1] + 2 * seqan.doc.constants.MARGIN;
    //console.log(t.dim);
    //console.log(seqan.doc.constants.MARGIN);

    // Draw hierarchy on canvas.
    seqan.doc.drawHierarchy(ctx, t);
    seqan.doc.createDivs('#class-links', t);
    $(window).resize(function() {
        seqan.doc.computeLayout(t, ctx);
        seqan.doc.updateDivs('#class-links', t);
    });
    //console.log(t)
}
