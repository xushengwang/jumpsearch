package Spiders::SpectrumViewer;

use strict;
use warnings;

use Carp;
use GD;
use GD::Graph::lines;
use GD::Graph::lines;
use GD::Graph::colour qw(:lists :colours);

our @ISA = qw(GD::Graph::Error);

our $VERSION = '0.01';


# Preloaded methods go here.

# Pairs of pattern and their associated annotations category used in
# conjunction with the colormap.  This is an array rather than a hash
# because we want to consider the patterns in a prescribed order.
my @defaultannotationsmatching = (
 '^y$|^y[ -]| y$| y[ -]', 'y',
 '^b$|^b[ -]| b$| b[ -]', 'b',
 '[iI]nternal', 'internal',
 '^$', 'none'
);


my %defaultcolormap = (
	'y' => 'red',
	'b' => 'blue',
	'internal' => 'green',
	'other' => 'dbrown',
	'none' => 'black'
);

my %fontmap = (
	'tiny' => gdTinyFont,
	'small' => gdSmallFont,
	'medium' => gdMediumBoldFont,
	'large' => gdLargeFont,
);

my %Defaults = (
	width => 500,
	height => 500,
	title => '',
	linewidth => 2,
	extranegativeheight => 0.1,
	ylabeldelta => 4, # offset of annotations on y axis measured in pixels
	xticknumber => 5, # of tick marks on X axis
	xlabeldelta => 6, # pixels for offset of annotations
	yaxismultiplier => 2.0, # a ratio, used to permit vertical room for peak annotations
	outputformat => 'png',
	peakfontsize => 'medium',
	x_label => 'm/z',
	y_label => 'Relative Intensity',
	annotationsmatching => \@defaultannotationsmatching,
	colormap => \%defaultcolormap,
);

sub new (\@\@;\@\%\%) # (masses, intensities, [optional] annotations, annotations_matching, colormap)
{
	my $type = shift;
	my $self = {};

	bless $self,$type;

	my($massesRef, $intensitiesRef, $annotationsRef, $annotations_matchingRef, $colormapRef) = @_;

	# Initialise all relevant parameters to defaults
	$self->initialise() or return;

	$self->{masses} = $massesRef;
	$self->{intensities} = $intensitiesRef;
	$self->{annotations} = $annotationsRef;
	$self->{annotationsmatching} = $annotations_matchingRef if $annotations_matchingRef;
	$self->{colormap} = $colormapRef if $colormapRef;

	return $self;
}

sub initialise
{
    my $self = shift;

    foreach (keys %Defaults) 
    {
        $self->set($_, $Defaults{$_});
    }

    return $self;
}


sub set {
        my ($self, $key, $value) = @_;

        $self->{$key} = $value;
}

sub plot
{
	my $self = shift;

	my $minmass;
	my $maxmass;
	my $maxintensity;
	my $minintensity;
	
	my @data_for_graph;
	
	my $i;
	my $j = 0;
	#
	# We are playing a trick by intentionally alternating the real y data values
	# with undefined values (and GD::Graph's skip_undef option), so that GD::Graph
	# doesn't try to plot our actual data ... we need to perform this plotting ourselves.
	# 
	for ($i = 0; $i <= $#{$self->{masses}}; $i++) {
		my $mass = $self->{masses}[$i];
		my $intensity = $self->{intensities}[$i];

		if ($i > 0) {
			$data_for_graph[0][$j] = $mass;
			$data_for_graph[1][$j] = undef;
			$j++;
		}
		$data_for_graph[0][$j] = $mass;
		$data_for_graph[1][$j] = $intensity;
		$minmass = $mass unless defined $minmass and $minmass < $mass;
		$maxmass = $mass unless defined $maxmass and $maxmass > $mass;
		$minintensity = $intensity unless defined $minintensity and $minintensity < $intensity;
		$maxintensity = $intensity unless defined $maxintensity and $maxintensity > $intensity;
		$j++;
	}
	
	#
	# adjust the min and max masses slightly, since otherwise
	# the min and max mass peaks will be obscured by the y axis and the graph's
	# right boundary
	#
	# we also force the ticks multiples of 5
	#
	my $massdiff = $maxmass - $minmass;
	my $tickconstraint = 5.0 * $self->{xticknumber};
	$minmass = int(($minmass - 0.04 * $massdiff)/$tickconstraint) * $tickconstraint;
	$maxmass = int(($maxmass + 0.04 * $massdiff)/$tickconstraint + 0.5) * $tickconstraint;
	$minmass = 0 if $minmass < 0;

	# note that we permit negative intensities; this permits some
	# interesting visualization capabilities
	#
	# extra vertical space is required to make the labels fit
	$minintensity = 0 if $minintensity > 0;
	$minintensity -= $self->{extranegativeheight} if $minintensity < 0;
	
	my $graph = GD::Graph::lines->new($self->{width},$self->{height});
	$graph->{graph}->setThickness($self->{linewidth});
	#
	# It turns out that if we don't specify x_tick_number explicitly, it becomes
	# quite messy to compute the conversion of x coordinates and requires the use
	# of lots of undocumented GD::Graph internals
	#
	# We are playing a trick with skip_undef, so that GD::Graph doesn't try to plot
	# our actual data ... we need to perform this plotting ourselves.
	# 
	# We double the maximum intensity so as to leave (hopefully) enough vertical
	# height for peak annotations
	#
	$graph->set(
		title               => $self->{title},
		x_label             => $self->{x_label},
		x_label_position    => 0.5,
		skip_undef          => 1,
		x_tick_number       => $self->{xticknumber},
		x_min_value         => $minmass,
		x_max_value         => $maxmass,
		x_number_format     => "%.1f",
		y_number_format     => "%.2f",
		y_min_value         => $minintensity * $self->{yaxismultiplier},
		y_max_value         => $maxintensity * $self->{yaxismultiplier},
		y_label => $self->{y_label}) or die $graph->error;
	
	$graph->set_x_axis_font(gdLargeFont);
	$graph->set_y_axis_font(gdLargeFont);

	#
	# draw the axes and their labels, and subsequently use the computed geometry
	# for scaling and translating our data points and annotations
	#
	my $im = $graph->plot(\@data_for_graph) or die $graph->error;
	
	my %colors;
	
	# make the background transparent and interlaced
	$colors{'white'} = $im->colorAllocate(255,255,255);
	$im->transparent($colors{'white'});
	$im->interlaced('true');
	
	for ($i = 0; $i <= $#{$self->{masses}}; $i++) {
		my $pattern;
		my $colorname;
		my $match = 'other';
		my $annot = $self->{annotations}[$i];
		my $mass = $self->{masses}[$i];
		my $intensity = $self->{intensities}[$i];
		my $discardThisAnnotation = 0;
		my $patternIndex;
	
		# by convention, a leading @ means discard this annotation, but
		# use it for purposes of coloring the peak
		if ($annot =~ m/^@/) {
			$discardThisAnnotation = 1;
			$annot =~ s/^.//;
		}
	
		PATTERN:
		for ($patternIndex = 0; $patternIndex < scalar(@{$self->{annotationsmatching}}); $patternIndex += 2) {
			$pattern = $self->{annotationsmatching}[$patternIndex];
			if (defined($annot) && $pattern && $annot =~ m/$pattern/) {
				$match = $self->{annotationsmatching}[$patternIndex + 1];
				last PATTERN;
			}
		}
		
		$colorname = $self->{colormap}{$match};
		_init_clr ($im, $colorname, \%colors) or carp "Unable to allocate color $colorname \n";
	#	print "match $match color $colorname \n";
		
		# draw vertical mass peaks and their annotations, if any
		_myline($self,$graph,$im,$mass,0,$mass,$intensity,$colors{$colorname});

		my $gdfont = $fontmap{$self->{peakfontsize}};
		$gdfont = gdMediumBoldFont unless $gdfont;

		# for negative values we label all mass peaks starting at the
		# bottom of the graph, since we lack the capability to
		# compute the vertical height of the labels and don't want
		# to require TrueType font availability in order to use
		# GD's stringFT method
		unless ($discardThisAnnotation) {
			if ($intensity >= 0) {
				_myannot($graph,$im,$mass,$intensity,$annot,$colors{$colorname},$self->{xlabeldelta},$self->{ylabeldelta},$gdfont);
			} else {
				_myannot($graph,$im,$mass,$minintensity*$self->{yaxismultiplier}*0.95,$annot,$colors{$colorname},$self->{xlabeldelta},0,$gdfont)
			}
		}
		
	}
	
	return $im->gif if ($self->{outputformat} eq 'gif');
	return $im->jpeg if ($self->{outputformat} eq 'jpeg');
	return $im->png;
	
}

sub _myline {
	my($self,$graph,$img,$xb,$yb,$xe,$ye,$color) = @_;
	my($xb2,$yb2,$xe2,$ye2);

	($xb2,$yb2) = $graph->val_to_pixel($xb,$yb,1);
	($xe2,$ye2) = $graph->val_to_pixel($xe,$ye,1);
	$img->line($xb2,$yb2,$xe2,$ye2,$color);
	$self->{_hotspots}{$xb} = ['line',$xb2,$yb2,$xe2,$ye2,$self->{linewidth}];
}

sub _myannot {
	my($graph,$img,$x,$y,$annot,$color,$xlabeldelta,$ylabeldelta,$font) = @_;
	my($x2,$y2) = $graph->val_to_pixel($x,$y,1);

	$img->stringUp($font,$x2-$xlabeldelta,$y2-$ylabeldelta,$annot,$color);
}

# swiped from set_clr in GD::Graph
sub _init_clr ($$\%) {
	my ($gd, $colorname, $colorsRef) = @_;

	my @rgb = _rgb($colorname);
	# All of this could potentially be done by using colorResolve
	# The problem is that colorResolve doesn't return an error
	# condition (-1) if it can't allocate a color. Instead it always
	# returns 0.
	
	# Check if this colour already exists on the canvas
	my $i = $gd->colorExact(@rgb);
	# if not, allocate a new one, and return its index
	$i = $gd->colorAllocate(@rgb) if $i < 0;
	# if this fails, we should use colorClosest.
	$i = $gd->colorClosest(@rgb)  if $i < 0;

	$$colorsRef{$colorname} = $i unless $i < 0;

	return $i;
	
}
1;
