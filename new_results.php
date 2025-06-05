<?php 
include 'header.php'; 
?>
<?php 
include 'db_connect.php'; 
?>
<main class="results-area">
    <h2>Analysis Results :</h2>
<?php
function clean_sequence($input) {
    $input = strtoupper($input);
    $input = preg_replace('/[^ATGCN]/', '', $input); 
    return $input;
}
function fetch_genbank_sequence($id) {
    $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$id&rettype=fasta&retmode=text";
    return @file_get_contents($url);
}
function fetch_uniprot_sequence($id) {
    $url = "https://rest.uniprot.org/uniprotkb/$id.fasta";
    return @file_get_contents($url);
}
$input = $_POST['sequence'] ?? '';
$input_type = $_POST['input_type'] ?? 'raw';
$sequence = '';
if ($input_type === 'fasta') {
    $lines = explode("\n", $input);
    foreach ($lines as $line) {
        $line = trim($line);
        if ($line && $line[0] !== '>') {
            $sequence .= $line;
        }
    }
} elseif ($input_type === 'genbank') {
    $fasta_data = fetch_genbank_sequence(trim($input));
    if ($fasta_data) {
        $lines = explode("\n", $fasta_data);
        foreach ($lines as $line) {
            if ($line && $line[0] !== '>') {
                $sequence .= trim($line);
            }
        }
    } else {
        echo "<p class='error'>Could not fetch sequence from GenBank.</p>";
    }
} elseif ($input_type === 'uniprot') {
    $fasta_data = fetch_uniprot_sequence(trim($input));
    if ($fasta_data) {
        $lines = explode("\n", $fasta_data);
        foreach ($lines as $line) {
            if ($line && $line[0] !== '>') {
                $sequence .= trim($line);
            }
        }
    } else {
        echo "<p class='error'>Could not fetch sequence from UniProt.</p>";
    }
} else {
    $sequence = $input;
}
$sequence = clean_sequence($sequence);
// echo $sequence;
$length = strlen($sequence);

$a = substr_count($sequence, 'A');
$t = substr_count($sequence, 'T');
$g = substr_count($sequence, 'G');
$c = substr_count($sequence, 'C');
$at = $a + $t;
$gc = $g + $c;
$at_percent = $length ? ($at / $length) * 100 : 0;
$gc_percent = $length ? ($gc / $length) * 100 : 0;

$molecular_weight = ($a * 313.21) + ($t * 304.2) + ($g * 329.21) + ($c * 289.18);
$melting_temp = $length ? (64.9 + 41 * ($gc - 16.4) / $length) : 0;

function complement($seq) {
    $complement = ['A' => 'T','T' => 'A','C' => 'G','G' => 'C', 'N' => 'N'];
    return implode('', array_map(fn($base) => $complement[$base] ?? 'N', str_split($seq)));
}
$complementary_strand = complement($sequence);
function iupac_to_regex($pattern) {
    $iupac = [
        'A' => 'A', 'C' => 'C', 'G' => 'G', 'T' => 'T',
        'R' => '[AG]', 'Y' => '[CT]', 'S' => '[GC]', 'W' => '[AT]',
        'K' => '[GT]', 'M' => '[AC]', 'B' => '[CGT]', 'D' => '[AGT]',
        'H' => '[ACT]', 'V' => '[ACG]', 'N' => '[ACGT]'
    ];
    return implode('', array_map(fn($b) => $iupac[$b] ?? $b, str_split(strtoupper($pattern))));
}
$matches = [];
$result = $conn->query("SELECT Enzyme_Name, Top_S, Bottom_S FROM y_6_pre_result");
while ($row = $result->fetch_assoc()) {
    $enzyme = $row['Enzyme_Name'];

    foreach (['Top_S' => 'top', 'Bottom_S' => 'bottom'] as $col => $strand) {
        $raw = $row[$col];
        // echo $raw;
        // echo "<h3>Enzyme: $enzyme ($strand strand)</h3>";
        $cut_offset = strpos($raw, '^');
        // echo "$cut_offset<br>";
        $site = str_replace('^', '', $raw);
        $regex = iupac_to_regex($site);
        $target_seq = ($strand === 'top') ? $sequence : $complementary_strand;

        if (preg_match_all("/$regex/", $target_seq, $all_matches, PREG_OFFSET_CAPTURE)) {
        foreach ($all_matches[0] as $match) {
        $start = $match[1];
        $cut_pos = $start + $cut_offset;

        if ($strand === 'top') {
            $cut_pos += 1;
        } else {
            $cut_pos = $length - $cut_pos+1; 
        }
// echo "$cut_pos<br>";
        $matches[] = [
           'enzyme' => $enzyme,
           'strand' => $strand,
           'position' => $cut_pos-1,
           'cut_seq' => $site,
           'match_start' => $start,
           'match_length' => strlen($site)
           ];
        }
     }
 }
}
?>
<div class="results-area">
  <div class="results-columns">
    <aside class="sequence-stats">
      <h2>Sequence Statistics</h2>
      <p><strong>Length:</strong> <?= $length ?> bp</p>
      <p><strong>A:</strong> <?= $a ?>  | <strong>T:</strong> <?= $t ?>  | <strong>G:</strong> <?= $g ?>  | <strong>C:</strong> <?= $c ?></p>
      <p><strong>AT Content:</strong> <?= round($at_percent, 2) ?>%</p>
      <p><strong>GC Content:</strong> <?= round($gc_percent, 2) ?>%</p>
      <p><strong>Molecular Weight:</strong> <?= round($molecular_weight, 2) ?> Da</p>
      <p><strong>Melting Temperature:</strong> <?= round($melting_temp, 2) ?> &deg;C</p>
    </aside>
  <main class="sequence-main">
      <h2>Input Sequence and Complement</h2>
      <pre style="font-family: monospace; white-space: pre;">
5'... <?= $sequence ?> ...3'
3'... <?= $complementary_strand ?> ...5'
</pre>
    </main>
<section class="cut-sites-section">
  <h2>Restriction Enzyme Cut Site Map</h2>
  <div style="margin-bottom:8px;">
    <button id="zoom-in" style="font-size:1.2em;">+</button>
    <button id="zoom-out" style="font-size:1.2em;">−</button>
    <span id="zoom-level" style="margin-left:10px;">100%</span>
  </div>
  <div id="svg-container" style="overflow-x:auto;">
    <svg id="enzyme-map" height="400" xmlns="http://www.w3.org/2000/svg" style="background: #fff; width: 100%; display: block;"></svg>
    <div id="enzyme-tooltip" style="display:none;position:absolute;pointer-events:none;background:#fff;border:1px solid #aaa;padding:8px 12px;font-size:13px;box-shadow:0 2px 8px #888;border-radius:6px;z-index:10"></div>
  </div>
</section>
<h2>Cut Site Table</h2>
<table border="1" cellpadding="8" cellspacing="0">
    <thead>
        <tr>
            <th>S. No.</th>
            <th>Enzyme Name</th>
            <th>Recognition Sequence</th>
            <th>Site Length</th>
            <th>Frequency</th>
            <th>Top Strand Cuts</th>
            <th>Bottom Strand Cuts</th>
        </tr>
    </thead>
    <tbody>
        <?php
        $enzyme_summary = [];
        foreach ($matches as $match) {
            $enzyme = $match['enzyme'];
            $cut_pos = $match['position'];
            $cut_seq = $match['cut_seq'];
        $total_rows = count($enzyme_summary);
            if (!isset($enzyme_summary[$enzyme])) {
               $enzyme_summary[$enzyme] = [
                    'sequence' => $cut_seq,
                    'length' => strlen($cut_seq),
                    'positions' => [],
                    'top' => [],
                    'bottom' => []
                ];
            }
            $enzyme_summary[$enzyme]['positions'][] = $cut_pos;
            $enzyme_summary[$enzyme][$match['strand']][] = $cut_pos;
        }
       $i = 1;
    foreach ($enzyme_summary as $enzyme => $info) {
    $rowClass = ($i > 25) ? " class='extra-row'" : "";
    echo "<tr{$rowClass}>";
    echo "<td>{$i}</td>";
    $enzyme_id = preg_replace('/\s+/', '_', $enzyme); // Clean ID
    echo "<td class='enzyme-cell' data-enzyme='{$enzyme_id}'>{$enzyme}</td>";
    echo "<td>{$info['sequence']}</td>";
    echo "<td>{$info['length']}</td>";
    echo "<td>" . count($info['positions']) . "</td>";
    // For top strand, subtract 1 from each position for display
    echo "<td>" . (empty($info['top']) ? '-' : implode(', ', array_map(fn($p) => $p - 1, $info['top']))) . "</td>";
    echo "<td>" . (empty($info['bottom']) ? '-' : implode(', ', $info['bottom'])) . "</td>";
    echo "</tr>";
    $i++;
}
        ?>
    </tbody>
</table>
<?php if ($total_rows > 25): ?>
    <button id="show-all-btn" class="show-all-btn">Show All</button>
<?php endif; ?>
<script>
const showBtn = document.getElementById("show-all-btn");
if (showBtn) {
  showBtn.addEventListener("click", function () {
    document.querySelectorAll(".extra-row").forEach(row => row.style.display = "table-row");
    this.style.display = "none";
  });
}
</script>
<script>
// --- 3. Visualize cut sites in a clean NEBcutter-like SVG map ---
const sequence = "<?= $sequence ?>";
const complement = "<?= complement($sequence) ?>";
const matches = <?= json_encode($matches) ?>;
const svg = document.getElementById("enzyme-map");
const tooltip = document.getElementById("enzyme-tooltip");
const charWidth = 18, charHeight = 20, padding = 30;
const seqY = 80, compY = 150;
svg.setAttribute("width", (charWidth * sequence.length + padding * 2));

// Draw sequence and ruler
function drawText(text, x, y, fontWeight = "normal") {
    const el = document.createElementNS("http://www.w3.org/2000/svg", "text");
    el.setAttribute("x", x); el.setAttribute("y", y);
    el.setAttribute("font-size", "15");
    el.setAttribute("font-family", "monospace");
    el.setAttribute("font-weight", fontWeight);
    el.textContent = text;
    svg.appendChild(el);
}
svg.innerHTML = ''; // clear
// Ruler ticks and numbers
for (let i = 1; i < sequence.length; i++) {
    const x = padding + i * charWidth;

    if (i % 10 === 0) {
        drawText(i.toString(), x+10, (seqY + compY) / 2, "normal");
        const tick = document.createElementNS("http://www.w3.org/2000/svg", "line");
        tick.setAttribute("x1", x+7);
        tick.setAttribute("y1", seqY+10); // Start just above top strand
        tick.setAttribute("x2", x+7);
        tick.setAttribute("y2", compY-20); // End just below bottom strand
        tick.setAttribute("stroke", "#888");
        tick.setAttribute("stroke-width", "0.1em");
        svg.appendChild(tick);
    }
    // Sequence letters
    drawText(sequence[i], x, seqY, "bold");
    drawText(complement[i], x, compY, "bold");
}
// Strand direction
drawText("5'... ", padding-35, seqY);
drawText("...3'", padding + sequence.length * charWidth + 5, seqY);
drawText("3'... ", padding-35, compY);
drawText("...5'", padding + sequence.length * charWidth + 5, compY);

// --- Group matches by cut position and strand for stacking ---
const grouped = {};
matches.forEach(m => {
    const key = m.strand + '_' + m.position;
    if (!grouped[key]) grouped[key] = [];
    grouped[key].push(m);
});

// --- Draw cut site lines and stacked enzyme labels ---
const labelBoxes = [];
Object.values(grouped).forEach(group => {
    group.sort((a, b) => a.enzyme.localeCompare(b.enzyme));
    const isTop = group[0].strand === 'top';
    // FIX for bottom strand: mirror cutX
    const cutX = isTop
        ? padding + group[0].position * charWidth
        : padding + (sequence.length - group[0].position) * charWidth;
    const baseY = isTop ? seqY-8 : compY+8;
    const stackBaseY = isTop ? seqY-35 : compY+35;
    const stackHeight = group.length * 16;
    const y2 = isTop ? (stackBaseY - stackHeight + 16) : (stackBaseY + stackHeight - 16);

    group.forEach((m, idx) => {
        let labelY = isTop ? (y2 + idx*16) : (y2 - (group.length-1-idx)*16 + 12);
        let hx = cutX + 12; // default: right side
        let placeLeft = false;

        // Check for overlap with previous label
        if (labelBoxes.length > 0) {
            const prev = labelBoxes[labelBoxes.length - 1];
            // If vertical distance is too small and horizontal overlap
            if (Math.abs(prev.y - labelY) < 16 && Math.abs(prev.x - (cutX + 12)) < 80) {
                // Check if there is blank space on left side
                if (cutX - 100 > 0) {
                    placeLeft = true;
                    // Shift label to left, but keep horizontal connector length the same as right side
                    // So, hx = cutX - (hx - cutX) = cutX - 12
                    hx = cutX - 12;
                } else {
                    // If not enough space left, and few enzymes, increase vertical arrow line length (move label further away)
                    if (group.length < 5) {
                        labelY += (idx+1) * 18; // stagger further down/up
                    }
                }
            }
        }

        // --- Draw vertical arrow line for this label ---
        const topArrowTip = seqY - 18;
        const topSeqY = seqY;
        const bottomSeqY = compY;
        const verticalSpace = topSeqY - topArrowTip; // typically 18
        const bottomArrowTip = bottomSeqY + verticalSpace;

        const arrowLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
        arrowLine.setAttribute("x1", cutX);
        arrowLine.setAttribute("x2", cutX);
        arrowLine.setAttribute("y1", isTop ? topArrowTip : bottomArrowTip);
        arrowLine.setAttribute("y2", labelY);
        arrowLine.setAttribute("stroke", "#e07a1f");
        arrowLine.setAttribute("stroke-width", "1");
        arrowLine.classList.add("enzyme-arrow-line");
        arrowLine.setAttribute("data-enzyme", m.enzyme.replace(/\s+/g, '_') + "_" + m.position + "_" + m.strand);
        svg.appendChild(arrowLine);

        // Arrow head (draw only for the first label in group to avoid duplicates)
        if (idx === 0) {
            const arrow = document.createElementNS("http://www.w3.org/2000/svg", "polygon");
            if (isTop) {
                arrow.setAttribute("points", `${cutX-6},${topArrowTip} ${cutX+6},${topArrowTip} ${cutX},${topArrowTip+10}`);
            } else {
                arrow.setAttribute("points", `${cutX-6},${bottomArrowTip} ${cutX+6},${bottomArrowTip} ${cutX},${bottomArrowTip-10}`);
            }
            arrow.setAttribute("fill", "#e07a1f");
            svg.appendChild(arrow);
        }

        // --- Draw horizontal connector ---
        const hLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
        hLine.setAttribute("x1", cutX);
        hLine.setAttribute("x2", hx);
        hLine.setAttribute("y1", labelY);
        hLine.setAttribute("y2", labelY);
        hLine.setAttribute("stroke", "#e07a1f");
        hLine.setAttribute("stroke-width", "1");
        hLine.classList.add("enzyme-hline");
        hLine.setAttribute("data-enzyme", m.enzyme.replace(/\s+/g, '_') + "_" + m.position + "_" + m.strand);
        svg.appendChild(hLine);

        // --- Enzyme label ---
        const text = document.createElementNS("http://www.w3.org/2000/svg", "text");
        text.setAttribute("x", placeLeft ? hx - 2 : hx + 2);
        text.setAttribute("y", labelY+4);
        text.setAttribute("font-size", "12");
        text.setAttribute("fill", "#e07a1f");
        text.setAttribute("font-family", "monospace");
        text.setAttribute("cursor", "pointer");
        text.classList.add("enzyme-label");
        const enzymeKey = m.enzyme.replace(/\s+/g, '_') + "_" + m.position + "_" + m.strand;
        text.setAttribute("data-enzyme", enzymeKey);
        text.setAttribute("data-cut", m.position);
        text.textContent = m.enzyme;
        svg.appendChild(text);

        // --- Tooltip and highlight logic ---
        text.addEventListener("mouseenter", function(e){
            highlightSeq(m, isTop, cutX);

            svg.querySelectorAll(`.enzyme-arrow-line[data-enzyme="${enzymeKey}"]`).forEach(line => {
                line.setAttribute("stroke", "#1976d2");
                line.setAttribute("stroke-width", "2.5");
            });
            svg.querySelectorAll(`.enzyme-hline[data-enzyme="${enzymeKey}"]`).forEach(line => {
                line.setAttribute("stroke", "#1976d2");
                line.setAttribute("stroke-width", "2.5");
            });

            tooltip.style.display = "block";
            const cutPos = isTop ? (m.position - 1) : (m.position );
            tooltip.innerHTML = `<b>${m.enzyme}</b> cuts at <b>${cutPos}</b><br>
                <b>Recognition sequence:</b> <span style="font-family:monospace;">${m.cut_seq}</span>`;
            const pt = svg.createSVGPoint();
            pt.x = Number(text.getAttribute("x")) + (text.getBBox().width || 0);
            pt.y = Number(text.getAttribute("y")) - 10;
            const ctm = text.getScreenCTM();
            if (ctm) {
                const transformed = pt.matrixTransform(ctm);
                tooltip.style.left = (window.scrollX + transformed.x + 10) + "px";
                tooltip.style.top = (window.scrollY + transformed.y) + "px";
            }
        });
        text.addEventListener("mouseleave", function(){
            tooltip.style.display = "none";
            removeHighlight();
            svg.querySelectorAll(".enzyme-arrow-line, .enzyme-hline").forEach(line => {
                line.setAttribute("stroke", "#e07a1f");
                line.setAttribute("stroke-width", "1");
            });
        });

        // Store for overlap checking
        labelBoxes.push({
            enzyme: m.enzyme,
            strand: m.strand,
            match_start: m.match_start,
            match_length: m.match_length,
            cut_pos: m.position,
            x: placeLeft ? hx - 2 : hx + 2,
            y: labelY+4,
            w: text.getBBox().width,
            h: 18,
            isTop: isTop
        });
    });
});

// --- Highlight sequence on hover ---
function highlightSeq(m, isTop, cutX) {
    removeHighlight();
    let xStart;
    if (isTop) {
        xStart = padding + m.match_start * charWidth;
    } else {
        // For bottom, highlight the actual matched region (not mirrored)
        xStart = padding + m.match_start * charWidth;
    }
    const y = m.strand === 'top' ? seqY-5 : compY-4;
    const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
    rect.setAttribute("x", xStart-4);
    rect.setAttribute("y", y-11);
    rect.setAttribute("width", m.match_length * charWidth);
    rect.setAttribute("height", 22);
    rect.setAttribute("fill", "#000000");
    rect.setAttribute("opacity", "0.5");
    rect.setAttribute("class", "highlight-box");
    // Insert highlight below text and labels
    svg.insertBefore(rect, svg.firstChild.nextSibling.nextSibling); // after background and before sequence text
}
function removeHighlight() {
    svg.querySelectorAll(".highlight-box").forEach(el => el.remove());
}

// --- Zoom in/out and pan ---
let zoom = 1;
let panX = 0, panY = 0;
const zoomInBtn = document.getElementById("zoom-in");
const zoomOutBtn = document.getElementById("zoom-out");
const zoomLevel = document.getElementById("zoom-level");
function setZoom() {
    svg.setAttribute("style", `background:#fff;width:100%;display:block;transform:scale(${zoom}) translate(${panX/zoom}px,${panY/zoom}px);transform-origin:0 0;`);
    zoomLevel.textContent = Math.round(zoom*100) + "%";
}
zoomInBtn.onclick = function() { zoom = Math.min(zoom*1.2, 5); setZoom(); };
zoomOutBtn.onclick = function() { zoom = Math.max(zoom/1.2, 0.2); setZoom(); };
svg.addEventListener('wheel', function(e){
    e.preventDefault();
    const oldZoom = zoom;
    if(e.deltaY < 0) zoom = Math.min(zoom*1.1, 5);
    else zoom = Math.max(zoom/1.1, 0.2);
    setZoom();
}, {passive:false});
// Pan with mouse drag
let isPanning = false, startX = 0, startY = 0, lastPanX = 0, lastPanY = 0;
svg.addEventListener('mousedown', function(e){
    isPanning = true;
    startX = e.clientX; startY = e.clientY;
    lastPanX = panX; lastPanY = panY;
});
window.addEventListener('mousemove', function(e){
    if(isPanning){
        panX = lastPanX + (e.clientX - startX);
        panY = lastPanY + (e.clientY - startY);
        setZoom();
    }
});
window.addEventListener('mouseup', function(){ isPanning = false; });
setZoom();
</script>

<?php include 'footer.php'; ?>