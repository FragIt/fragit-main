$LOADCOMMAND

function ColorDefault() {
  select all
  color cpk
  select none
}

function ColorFragments(type) {
  if type == "default" {
    ColorDefault()
  }

  if type == "fragments" {
    ColorByFragment()
  }
  if type == "buffer" {
    ColorByBuffer()
  }
  if type == "active" {
    ColorByActive()
  }
}
function ColorByFragment() {
$FRAGMENTS
select none
}
function ColorByBuffer() {
select all
color green
select none
$BUFFER
select none
}
function ColorByActive() {
select all
color green
select none
$BUFFER
$ACTIVE
select none
}

function ColorByBackbone() {
  ColorDefault()
$BACKBONE
select none
}

function LabelFragments() {
select none
}

function DrawBreakPoints() {
$BREAKPOINTS
}

function HideBreakPoints() {
draw cb* off
}

function ShowBreakPoints() {
draw cb* on
}

ColorFragments("fragments")
DrawBreakPoints()
HideBreakPoints()
