# Please developers (not translators) don't reuse a key inside another key
# or concat text with a translation string in the code, use the appropriate
# features like variables and selectors in each appropriate case!
# And also comment the string if it's something not easy to translate, so you help
# translators on what it means


## Websocket (server) status

websocket-connecting = Verbinding maken met de server
websocket-connection_lost = Verbinding met de server verbroken. Opniew verbinding maken...
websocket-connection_lost-desc = Het ziet er naar uit dat de verbinding met de SlimeVR server is verbroken. Check het logboek en start het programma opnieuw.
websocket-timedout = Kan niet verbinden met de server.
websocket-timedout-desc = Het ziet er naar uit dat de SlimeVR server is gestopt. Check het logboek en start het programma opnieuw.
websocket-error-close = SlimeVR afsluiten
websocket-error-logs = Open het logboek.

## Update notification

version_update-title = Nieuwe versie beschikbaar: { $version }
version_update-description = Als je op "{ version_update-update }" klikt, wordt het SlimeVR-installatieprogramma voor je gedownload.
version_update-update = Bijwerken
version_update-close = Sluiten

## Tips

tips-find_tracker = Weet je niet welke tracker welke is? Schud een tracker en het corresponderende item zal worden gemarkeerd.
tips-do_not_move_heels = Zorg ervoor dat je hielen niet bewegen tijdens de opname!
tips-file_select = Sleep bestanden naar hier om ze te gebruiken of <u>blader</u>.
tips-tap_setup = Je kan langzaam 2 keer op je tracker tikken om deze te kiezen in plaats van deze in het menu te selecteren.
tips-turn_on_tracker = Gebruik je officiële SlimeVR-trackers? Vergeet niet om <b><em>je tracker aan te zetten</em></b> nadat je deze op de pc hebt aangesloten!
tips-failed_webgl = WebGL initialiseren is gefaald.

## Body parts

body_part-NONE = Niet toegewezen
body_part-HEAD = Hoofd
body_part-NECK = Nek
body_part-RIGHT_SHOULDER = Rechterschouder
body_part-RIGHT_UPPER_ARM = Rechterbovenarm
body_part-RIGHT_LOWER_ARM = Rechteronderarm
body_part-RIGHT_HAND = Rechterhand
body_part-RIGHT_UPPER_LEG = Rechterdij
body_part-RIGHT_LOWER_LEG = Rechterenkel
body_part-RIGHT_FOOT = Rechtervoet
body_part-UPPER_CHEST = Bovenborst
body_part-CHEST = Borst
body_part-WAIST = Taille
body_part-HIP = Heup
body_part-LEFT_SHOULDER = Linkerschouder
body_part-LEFT_UPPER_ARM = Linkerbovenarm
body_part-LEFT_LOWER_ARM = Linkeronderarm
body_part-LEFT_HAND = Linkerhand
body_part-LEFT_UPPER_LEG = Linkerdij
body_part-LEFT_LOWER_LEG = Linkerenkel
body_part-LEFT_FOOT = Linkervoet
body_part-LEFT_THUMB_METACARPAL = Linkerduim middenhandsbeentje
body_part-LEFT_THUMB_PROXIMAL = Linkerduim proximaal
body_part-LEFT_THUMB_DISTAL = Linkerduim distaal
body_part-LEFT_INDEX_PROXIMAL = Linker wijsvinger proximaal
body_part-LEFT_INDEX_INTERMEDIATE = Linker middelste kootje van de wijsvinger
body_part-LEFT_INDEX_DISTAL = Linker wijsvinger distaal
body_part-LEFT_MIDDLE_PROXIMAL = Linker middelvinger proximaal
body_part-LEFT_MIDDLE_INTERMEDIATE = Linker middelste kootje van de middelvinger
body_part-LEFT_MIDDLE_DISTAL = Linker middelvinger distaal
body_part-LEFT_RING_PROXIMAL = Linker ringvinger proximaal
body_part-LEFT_RING_INTERMEDIATE = Linker middelste kootje van de ringvinger
body_part-LEFT_RING_DISTAL = Linker ringvinger distaal
body_part-LEFT_LITTLE_PROXIMAL = Linker kleine vinger proximaal
body_part-LEFT_LITTLE_INTERMEDIATE = Linker middelste kootje van de kleine vinger
body_part-LEFT_LITTLE_DISTAL = Linker kleine vinger distaal
body_part-RIGHT_THUMB_METACARPAL = Rechterduim middenhandsbeentje
body_part-RIGHT_THUMB_PROXIMAL = Rechterduim proximaal
body_part-RIGHT_THUMB_DISTAL = Rechterduim distaal
body_part-RIGHT_INDEX_PROXIMAL = Rechter wijsvinger proximaal
body_part-RIGHT_INDEX_INTERMEDIATE = Rechter middelste kootje van de wijsvinger
body_part-RIGHT_INDEX_DISTAL = Rechter wijsvinger distaal
body_part-RIGHT_MIDDLE_PROXIMAL = Rechts middelvinger proximaal
body_part-RIGHT_MIDDLE_INTERMEDIATE = Rechter middelste kootje van de middelvinger
body_part-RIGHT_MIDDLE_DISTAL = Rechter middelvinger distaal
body_part-RIGHT_RING_PROXIMAL = Rechter ringvinger proximaal
body_part-RIGHT_RING_INTERMEDIATE = Rechter middelste kootje van de ringvinger
body_part-RIGHT_RING_DISTAL = Rechter ringvinger distaal
body_part-RIGHT_LITTLE_PROXIMAL = Rechter kleine vinger proximaal
body_part-RIGHT_LITTLE_INTERMEDIATE = Rechter middelste kootje van de kleine vinger
body_part-RIGHT_LITTLE_DISTAL = Rechter kleine vinger distaal

## BoardType

board_type-UNKNOWN = Onbekend
board_type-NODEMCU = NodeMCU
board_type-CUSTOM = Custom bord
board_type-WROOM32 = WROOM32
board_type-WEMOSD1MINI = Wemos D1 Mini
board_type-TTGO_TBASE = TTGO T-Base
board_type-ESP01 = ESP-01
board_type-SLIMEVR = SlimeVR
board_type-LOLIN_C3_MINI = Lolin C3 Mini
board_type-BEETLE32C3 = Beetle ESP32-C3
board_type-ESP32C3DEVKITM1 = Espressif ESP32-C3 DevKitM-1
board_type-OWOTRACK = owoTrack
board_type-WRANGLER = Wrangler Joycons
board_type-MOCOPI = Sony Mocopi
board_type-WEMOSWROOM02 = Wemos Wroom-02 D1 Mini
board_type-XIAO_ESP32C3 = Seeed Studio XIAO ESP32C3
board_type-HARITORA = Haritora
board_type-ESP32C6DEVKITC1 = Espressif ESP32-C6 DevKitC-1
board_type-GLOVE_IMU_SLIMEVR_DEV = SlimeVR Dev IMU Handschoen

## Proportions

skeleton_bone-NONE = Geen
skeleton_bone-HEAD = Hoofdverschuiving
skeleton_bone-HEAD-desc =
    Dit is de afstand tussen je headset en het midden van je hoofd.
    Om dit aan te passen, schud je je hoofd naar links en rechts alsof je 'nee' zegt, 
    pas het aan totdat beweging van de andere trackers te verwaarlozen is.
skeleton_bone-NECK = Neklengte
skeleton_bone-NECK-desc =
    Dit is de afstand tussen het midden van je hoofd en de basis van je nek.
    Om dit aan te passen, beweeg je je hoofd op en neer alsof je knikt, of kantel je je hoofd 
    naar links en rechts. Wijzig de positie totdat beweging in andere trackers verwaarloosbaar is.
skeleton_bone-torso_group = Romp lengte
skeleton_bone-torso_group-desc =
    Dit is de afstand van je nek tot je heupen.
    Om dit aan te passen, ga rechtop staan en pas het aan totdat je virtuele heupen
    in lijn zijn met je echte heupen.
skeleton_bone-UPPER_CHEST = Bovenborst Lengte
skeleton_bone-UPPER_CHEST-desc =
    Dit is de afstand tussen de basis van je nek en het midden van je borst.
    Om dit aan te passen, stel je de torso-lengte correct af en pas je deze aan in verschillende houdingen 
    (zitten, bukken, liggen, enz.) totdat je virtuele ruggengraat overeenkomt met je echte.
skeleton_bone-CHEST_OFFSET = Borstoffset
skeleton_bone-CHEST_OFFSET-desc =
    Dit kan worden aangepast om je virtuele borsttracker omhoog of omlaag te verplaatsen,
    om te helpen bij de kalibratie in bepaalde spellen of applicaties die verwachten dat deze hoger of lager staat.
skeleton_bone-CHEST = Borstafstand
skeleton_bone-CHEST-desc =
    Dit is de afstand van het midden van je borst tot het midden van je ruggengraat.
    Om dit aan te passen, stel je de torso-lengte correct af en pas je deze aan in verschillende houdingen 
    (zitten, bukken, liggen, enz.) totdat je virtuele ruggengraat overeenkomt met je echte.
skeleton_bone-WAIST = Taille lengte
skeleton_bone-WAIST-desc =
    Dit is de afstand van het midden van je ruggengraat tot je navel.
    Om dit aan te passen, stel je de torso-lengte correct af en pas je deze aan in verschillende houdingen 
    (zitten, bukken, liggen, enz.) totdat je virtuele ruggengraat overeenkomt met je echte.
skeleton_bone-HIP = Heuplengte
skeleton_bone-HIP-desc =
    Dit is de afstand van je navel tot je heupen.
    Om dit aan te passen, stel je de torso-lengte correct in en pas je deze aan in verschillende houdingen
    (zitten, bukken, liggen, enz.) totdat je virtuele ruggengraat overeenkomt met je echte.
skeleton_bone-HIP_OFFSET = Heupoffset
skeleton_bone-HIP_OFFSET-desc =
    Dit kan worden aangepast om je virtuele heuptracker omhoog of omlaag te verplaatsen,
    om te helpen bij de kalibratie in bepaalde spellen of applicaties die mogelijk verwachten dat deze zich rond je middel bevindt.
skeleton_bone-HIPS_WIDTH = Heupbreedte
skeleton_bone-HIPS_WIDTH-desc =
    Dit is de afstand tussen het begin van je benen.
    Om dit aan te passen, voer je een volledige reset uit met je benen gestrekt en pas je het aan totdat je virtuele benen horizontaal overeenkomen met je echte.
skeleton_bone-leg_group = Beenlengte
skeleton_bone-leg_group-desc =
    Dit is de afstand van je heupen tot je voeten.
    Om dit aan te passen, pas je je torso-lengte op de juiste manier aan
    totdat je virtuele voeten op hetzelfde niveau staan als je echte.
skeleton_bone-UPPER_LEG = Bovenbeenlengte
skeleton_bone-UPPER_LEG-desc =
    Dit is de afstand van je heupen tot je knieën.
    Om dit aan te passen, pas je je beenlengte op de juiste manier aan
    totdat je virtuele knieën op dezelfde hoogte zijn als je echte.
skeleton_bone-LOWER_LEG = Onderbeenlengte
skeleton_bone-LOWER_LEG-desc =
    Dit is de afstand van je knieën tot je enkels.
    Om dit aan te passen, pas je je beenlengte op de juiste manier aan
    totdat je virtuele knieën op dezelfde hoogte zijn als je echte knieën.
skeleton_bone-FOOT_LENGTH = Voetlengte
skeleton_bone-FOOT_LENGTH-desc =
    Dit is de afstand van je enkels tot je tenen.
    Om dit aan te passen, ga op je tenen staan en pas het aan totdat je virtuele voeten op hun plaats blijven.
skeleton_bone-FOOT_SHIFT = Voetverschuiving
skeleton_bone-FOOT_SHIFT-desc =
    Deze waarde is de horizontale afstand van je knie tot je enkel.
    Dit houdt rekening met het feit dat je onderbenen naar achteren staan wanneer je rechtop staat.
    Om dit aan te passen, stel je de voetlengte in op 0, voer je een volledige reset uit,
    en pas je het aan totdat je virtuele voeten op één lijn liggen met het midden van je enkels.
skeleton_bone-SKELETON_OFFSET = Skelet offset
skeleton_bone-SKELETON_OFFSET-desc =
    Dit kan worden aangepast om al je trackers naar voren of naar achteren te verschuiven.
    Het kan worden gebruikt om te helpen bij de kalibratie in bepaalde spellen of toepassingen
    die mogelijk verwachten dat je trackers verder naar voren staan.
skeleton_bone-SHOULDERS_DISTANCE = Schoudersafstand
skeleton_bone-SHOULDERS_DISTANCE-desc =
    Dit is de verticale afstand van de basis van je nek tot je schouders.
    Om dit aan te passen, stel je de lengte van je bovenarm in op 0 en
    pas je deze aan totdat je virtuele elleboogtrackers verticaal uitlijnen met je echte schouders.
skeleton_bone-SHOULDERS_WIDTH = Schouderbreedte
skeleton_bone-SHOULDERS_WIDTH-desc =
    Dit is de horizontale afstand van de basis van je nek tot je schouders.
    Om dit aan te passen, stel je de lengte van je bovenarm in op 0 en pas je deze aan
    totdat je virtuele elleboogtrackers horizontaal uitlijnen met je echte schouders.
skeleton_bone-arm_group = Armlengte
skeleton_bone-arm_group-desc =
    Dit is de afstand van je schouders tot je polsen.
    Om dit aan te passen, pas je de schouderafstand correct aan, stel je Handafstand Y in op 0, 
    en pas je deze aan totdat je handtrackers op één lijn liggen met je polsen.
skeleton_bone-UPPER_ARM = Bovenarmlengte
skeleton_bone-UPPER_ARM-desc =
    Dit is de afstand van je schouders tot je ellebogen.
    Om dit aan te passen, pas je de armlengte correct aan en pas je deze aan
    totdat je elleboogtrackers overeenkomen met je echte ellebogen.
skeleton_bone-LOWER_ARM = Onderarmlengte
skeleton_bone-LOWER_ARM-desc =
    Dit is de afstand van je ellebogen tot je polsen.
    Om dit aan te passen, pas je de armlengte correct aan en pas je deze aan
    totdat je elleboogtrackers overeenkomen met je echte ellebogen.
skeleton_bone-HAND_Y = Afstand hand Y
skeleton_bone-HAND_Y-desc =
    Dit is de verticale afstand van je polsen tot het midden van je hand.
    Om dit aan te passen voor motion capture, pas je de armlengte correct aan 
    en pas je deze aan totdat je handtrackers verticaal uitgelijnd zijn met het midden van je handen. 
    Wil je het aanpassen voor elleboogtracking vanaf je controllers, 
    stel dan de armlengte in op 0 en pas je deze aan totdat je elleboogtrackers verticaal op één lijn liggen met je polsen.
skeleton_bone-HAND_Z = Afstand hand Z
skeleton_bone-HAND_Z-desc =
    Dit is de horizontale afstand van je polsen tot het midden van je hand.
    Als je dit wilt aanpassen voor motion capture, stel je deze in op 0.
    Wil je het aanpassen voor elleboogtracking vanaf je controllers, stel dan de armlengte in op 0 
    en pas je deze aan totdat je elleboogtrackers horizontaal op één lijn liggen met je polsen.
skeleton_bone-ELBOW_OFFSET = Elleboogoffset
skeleton_bone-ELBOW_OFFSET-desc = Dit kan worden aangepast om je virtuele elleboogtrackers omhoog of omlaag te verplaatsen, zodat wordt voorkomen dat VRChat per ongeluk een elleboogtracker aan de borst koppelt.

## Tracker reset buttons

reset-reset_all = Alle afmetingen resetten
reset-reset_all_warning-v2 =
    <b>Waarschuwing:</b> Je verhoudingen worden teruggezet naar de standaardinstellingen, geschaald op basis van je ingestelde lengte.
    Weet je zeker dat je dit wilt doen?
reset-reset_all_warning-reset = Verhoudingen resetten
reset-reset_all_warning-cancel = Annuleren
reset-reset_all_warning_default-v2 =
    <b>Waarschuwing:</b> Uw lengte is nog niet ingesteld, je verhoudingen worden teruggezet naar de standaardinstellingen met de standaard lengte.
    Weet je zeker dat je dit wilt doen?
reset-full = Volledige reset
reset-mounting = Reset montage
reset-mounting-feet = Reset voetmontage
reset-mounting-fingers = Reset vingermontage
reset-yaw = Yaw Reset

## Serial detection stuff

serial_detection-new_device-p0 = Nieuw serieel apparaat gedetecteerd!
serial_detection-new_device-p1 = Voer je WiFi-inloggegevens in!
serial_detection-new_device-p2 = Selecteer wat je wil doen
serial_detection-open_wifi = Verbind met Wi-Fi
serial_detection-open_serial = Seriële console openen
serial_detection-submit = Verzenden!
serial_detection-close = Sluiten

## Navigation bar

navbar-home = Startpagina
navbar-body_proportions = Lichaamsverhoudingen
navbar-trackers_assign = Tracker-toewijzing
navbar-mounting = Montage-kalibratie
navbar-onboarding = Installatiewizard
navbar-settings = Instellingen

## Biovision hierarchy recording

bvh-start_recording = BVH opnemen
bvh-recording = Opname bezig...
bvh-save_title = Sla BVH-opname op

## Tracking pause

tracking-unpaused = Pauzeer tracking
tracking-paused = Hervat tracking

## Widget: Overlay settings

widget-overlay = Overlay
widget-overlay-is_visible_label = Overlay in SteamVR weergeven
widget-overlay-is_mirrored_label = Overlay weergeven als spiegel

## Widget: Drift compensation

widget-drift_compensation-clear = Reset huidige drift compensatie

## Widget: Clear Reset Mounting

widget-clear_mounting = Reset montage legen

## Widget: Developer settings

widget-developer_mode = Developer Mode
widget-developer_mode-high_contrast = Hoog contrast
widget-developer_mode-precise_rotation = Precieze rotatie
widget-developer_mode-fast_data_feed = Snelle data feed
widget-developer_mode-filter_slimes_and_hmd = Filter slimes en HMD
widget-developer_mode-sort_by_name = Sorteer op naam
widget-developer_mode-raw_slime_rotation = Ruwe rotatie
widget-developer_mode-more_info = Meer informatie

## Widget: IMU Visualizer

widget-imu_visualizer = Rotatie
widget-imu_visualizer-preview = Voorbeeld
widget-imu_visualizer-hide = Verbergen
widget-imu_visualizer-rotation_raw = Rauw
widget-imu_visualizer-rotation_preview = Preview
widget-imu_visualizer-acceleration = Versnelling
widget-imu_visualizer-position = Positie
widget-imu_visualizer-stay_aligned = Blijf in lijn

## Widget: Skeleton Visualizer

widget-skeleton_visualizer-preview = Skelet voorbeeld
widget-skeleton_visualizer-hide = Verbergen

## Tracker status

tracker-status-none = Geen status
tracker-status-busy = Bezig
tracker-status-error = Fout
tracker-status-disconnected = Verbinding verbroken
tracker-status-occluded = Verborgen
tracker-status-ok = OK
tracker-status-timed_out = Timed Out

## Tracker status columns

tracker-table-column-name = Naam
tracker-table-column-type = Type
tracker-table-column-battery = Batterij
tracker-table-column-ping = Ping
tracker-table-column-tps = TPS
tracker-table-column-temperature = Temp. °C
tracker-table-column-linear-acceleration = Accel. X/Y/Z
tracker-table-column-rotation = Rotatie X/Y/Z
tracker-table-column-position = Positie X/Y/Z
tracker-table-column-stay_aligned = Blijf in lijn
tracker-table-column-url = URL

## Tracker rotation

tracker-rotation-front = Voorzijde
tracker-rotation-front_left = Linksvoor
tracker-rotation-front_right = Rechtsvoor
tracker-rotation-left = Links
tracker-rotation-right = Rechts
tracker-rotation-back = Achterzijde
tracker-rotation-back_left = Linksachter
tracker-rotation-back_right = Rechtsachter
tracker-rotation-custom = Aangepast
tracker-rotation-overriden = (overschreven door montage reset)

## Tracker information

tracker-infos-manufacturer = Fabrikant
tracker-infos-display_name = Weergavenaam
tracker-infos-custom_name = Aangepaste naam
tracker-infos-url = Tracker URL
tracker-infos-version = Firmware versie
tracker-infos-hardware_rev = Hardware revisie
tracker-infos-hardware_identifier = Hardware-id
tracker-infos-data_support = Gegevensondersteuning
tracker-infos-imu = IMU-sensor
tracker-infos-board_type = Mainbord
tracker-infos-network_version = Protocol Versie
tracker-infos-magnetometer = Magnetometer
tracker-infos-magnetometer-status-v1 =
    { $status ->
        [DISABLED] Uitgeschakeld
        [ENABLED] Ingeschakeld
       *[NOT_SUPPORTED] Niet ondersteund
    }

## Tracker settings

tracker-settings-back = Terug naar trackerslijst
tracker-settings-title = Trackersinstellingen
tracker-settings-assignment_section = Toewijzing
tracker-settings-assignment_section-description = Aan welk lichaamsdeel de tracker is toegewezen.
tracker-settings-assignment_section-edit = Toewijzing bewerken
tracker-settings-mounting_section = Montage oriëntatie
tracker-settings-mounting_section-description = Waar is de tracker gemonteerd?
tracker-settings-mounting_section-edit = Montage bewerken
tracker-settings-drift_compensation_section = Laat drift compensatie toe
tracker-settings-drift_compensation_section-description = Moet deze tracker compenseren voor drift wanneer drift compensatie is ingeschakeld?
tracker-settings-drift_compensation_section-edit = Laat drift compensatie toe
tracker-settings-use_mag = Sta de magnetometer toe op deze tracker.
# Multiline!
tracker-settings-use_mag-description =
    Wilt u dat deze tracker de magnetometer gebruikt om drift te verminderen wanneer de magnetometer is toegestaan? <b>Zet de tracker niet uit terwijl u dit aan of uit zet.</b>
    U moet eerst de magnetometer toestemming geven,<magSetting>click hier om naar de instellingen te gaan</magSetting>.
tracker-settings-use_mag-label = Laat magnetometer toe
# The .<name> means it's an attribute and it's related to the top key.
# In this case that is the settings for the assignment section.
tracker-settings-name_section = Trackernaam
tracker-settings-name_section-description = Geef een schattige bijnaam :)
tracker-settings-name_section-placeholder = NightyBeast's linkerbeen
tracker-settings-name_section-label = Trackernaam
tracker-settings-forget = Vergeet tracker
tracker-settings-forget-description = Verwijdert de tracker van de SlimeVR Server en voorkomt dat deze verbinding kan maken totdat de server opnieuw wordt opgestart. De configuratie van de tracker blijft behouden.
tracker-settings-forget-label = Vergeet tracker
tracker-settings-update-unavailable = Kan niet worden bijgewerkt (DIY)
tracker-settings-update-low-battery = Kan niet worden bijgewerkt. Batterij lager dan 50%
tracker-settings-update-up_to_date = Up to date.
tracker-settings-update-blocked = Update is niet beschikbaar. Er zijn geen andere versies beschikbaar.
tracker-settings-update-available = { $versionName } is nu beschikbaar
tracker-settings-update = Werk nu bij.
tracker-settings-update-title = Firmware versie

## Tracker part card info

tracker-part_card-no_name = Geen naam
tracker-part_card-unassigned = Niet toegewezen

## Body assignment menu

body_assignment_menu = Waar wil je deze tracker bevestigen?
body_assignment_menu-description = Kies een locatie waar je deze tracker wilt toewijzen. Als alternatief kun je kiezen om alle trackers tegelijk te beheren in plaats van één voor één.
body_assignment_menu-show_advanced_locations = Geavanceerde bevestigingslocaties weergeven
body_assignment_menu-manage_trackers = Beheer alle trackers
body_assignment_menu-unassign_tracker = Tracker niet toewijzen

## Tracker assignment menu

# A -translation_key (with a dash in the front) means that it's a label.
# It can only be used in the translation file, it's nice for reusing names and that kind of stuff.
#
# We are using it here because english doesn't require changing the text in each case but
# maybe your language does.
-tracker_selection-part = Welke tracker wil je toewijzen aan je
tracker_selection_menu-NONE = Van welke tracker will je de toewijzing ongedaan maken?
tracker_selection_menu-HEAD = { -tracker_selection-part } hoofd?
tracker_selection_menu-NECK = { -tracker_selection-part } nek?
tracker_selection_menu-RIGHT_SHOULDER = { -tracker_selection-part } rechterschouder?
tracker_selection_menu-RIGHT_UPPER_ARM = { -tracker_selection-part } rechterbovenarm?
tracker_selection_menu-RIGHT_LOWER_ARM = { -tracker_selection-part } rechteronderarm?
tracker_selection_menu-RIGHT_HAND = { -tracker_selection-part } rechterhand?
tracker_selection_menu-RIGHT_UPPER_LEG = { -tracker_selection-part } rechterdij?
tracker_selection_menu-RIGHT_LOWER_LEG = { -tracker_selection-part } rechterenkel?
tracker_selection_menu-RIGHT_FOOT = { -tracker_selection-part } rechtervoet?
tracker_selection_menu-RIGHT_CONTROLLER = { -tracker_selection-part } rechtercontroller?
tracker_selection_menu-UPPER_CHEST = { -tracker_selection-part } bovenborst?
tracker_selection_menu-CHEST = { -tracker_selection-part } borst?
tracker_selection_menu-WAIST = { -tracker_selection-part } taille?
tracker_selection_menu-HIP = { -tracker_selection-part } heup?
tracker_selection_menu-LEFT_SHOULDER = { -tracker_selection-part } linkerschouder?
tracker_selection_menu-LEFT_UPPER_ARM = { -tracker_selection-part } linkerbovenarm?
tracker_selection_menu-LEFT_LOWER_ARM = { -tracker_selection-part } linkeronderarm?
tracker_selection_menu-LEFT_HAND = { -tracker_selection-part } linkerhand?
tracker_selection_menu-LEFT_UPPER_LEG = { -tracker_selection-part } linkerdij?
tracker_selection_menu-LEFT_LOWER_LEG = { -tracker_selection-part } linkerenkel?
tracker_selection_menu-LEFT_FOOT = { -tracker_selection-part } linkervoet?
tracker_selection_menu-LEFT_CONTROLLER = { -tracker_selection-part } linkercontroller?
tracker_selection_menu-unassigned = Niet toegewezen trackers
tracker_selection_menu-assigned = Toegewezen trackers
tracker_selection_menu-dont_assign = Niet toewijzen
# This line cares about multilines.
# <b>text</b> means that the text should be bold.
tracker_selection_menu-neck_warning =
    <b>Waarschuwing:</b> Een nektracker kan dodelijk zijn indien deze te strak wordt afgesteld,
    de band kan de bloedsomloop naar je hoofd afsnijden!
tracker_selection_menu-neck_warning-done = Ik begrijp de risico's
tracker_selection_menu-neck_warning-cancel = Annuleren

## Mounting menu

mounting_selection_menu = Waar wil je deze tracker hebben bevestigd?
mounting_selection_menu-close = Sluiten

## Sidebar settings

settings-sidebar-title = Instellingen
settings-sidebar-general = Algemeen
settings-sidebar-tracker_mechanics = Trackersinstellingen
settings-sidebar-stay_aligned = Blijf in lijn
settings-sidebar-fk_settings = FK-instellingen
settings-sidebar-gesture_control = Tikbediening
settings-sidebar-interface = Interface
settings-sidebar-osc_router = OSC-router
settings-sidebar-osc_trackers = VRChat OSC Trackers
settings-sidebar-utils = Hulpmiddelen
settings-sidebar-serial = Serieel console
settings-sidebar-appearance = Uiterlijk
settings-sidebar-notifications = Notificaties
settings-sidebar-behavior = Gedrag
settings-sidebar-firmware-tool = DIY Firmware Tool
settings-sidebar-vrc_warnings = VRChat Configuratie-waarschuwingen
settings-sidebar-advanced = Geavanceerd

## SteamVR settings

settings-general-steamvr = SteamVR
settings-general-steamvr-subtitle = SteamVR trackers
# Not all translation keys support multiline, only the ones that specify it will actually
# split it in lines (that also means you can split in lines however you want in those).
# The first spaces (not tabs) for indentation will be ignored, just to make the file look nice when writing.
# This one is one of this cases that cares about multilines
settings-general-steamvr-description =
    Schakel specifieke SteamVR trackers in of uit.
    Handig voor games of apps die alleen bepaalde trackers ondersteunen.
settings-general-steamvr-trackers-waist = Taille
settings-general-steamvr-trackers-chest = Borst
settings-general-steamvr-trackers-left_foot = Linkervoet
settings-general-steamvr-trackers-right_foot = Rechtervoet
settings-general-steamvr-trackers-left_knee = Linkerknie
settings-general-steamvr-trackers-right_knee = Rechterknie
settings-general-steamvr-trackers-left_elbow = Linker elleboog
settings-general-steamvr-trackers-right_elbow = Rechter elleboog
settings-general-steamvr-trackers-left_hand = Linkerhand
settings-general-steamvr-trackers-right_hand = Rechterhand
settings-general-steamvr-trackers-tracker_toggling = Automatische tracker toewijzing
settings-general-steamvr-trackers-tracker_toggling-description = Zorgt automatisch voor het in- en uitschakelen van SteamVR-trackers, afhankelijk van je huidige tracker toewijzingen.
settings-general-steamvr-trackers-tracker_toggling-label = Automatische tracker toewijzing
settings-general-steamvr-trackers-hands-warning = <b>Waarschuwing:</b> hand trackers negeren je controllers. Weet je zeker dat je wilt doorgaan?
settings-general-steamvr-trackers-hands-warning-cancel = Annuleren
settings-general-steamvr-trackers-hands-warning-done = Ja

## Tracker mechanics

settings-general-tracker_mechanics = Tracker aanpassingen
settings-general-tracker_mechanics-filtering = Filtering
# This also cares about multilines
settings-general-tracker_mechanics-filtering-description =
    Kies het type filter voor je trackers.
    Voorspelling voorspelt beweging terwijl smoothing bewegingen vloeiender maakt.
settings-general-tracker_mechanics-filtering-type = Filtering type
settings-general-tracker_mechanics-filtering-type-none = Geen filtering
settings-general-tracker_mechanics-filtering-type-none-description = Gebruik rotaties zoals ze zijn. Zal geen filtering uitvoeren.
settings-general-tracker_mechanics-filtering-type-smoothing = Smoothing
settings-general-tracker_mechanics-filtering-type-smoothing-description = Maakt bewegingen vloeiender, maar voegt enige latentie toe.
settings-general-tracker_mechanics-filtering-type-prediction = Voorspelling
settings-general-tracker_mechanics-filtering-type-prediction-description = Verlaagt latentie en maakt bewegingen snappier, maar kan jitter verhogen.
settings-general-tracker_mechanics-filtering-amount = Hoeveelheid
settings-general-tracker_mechanics-yaw-reset-smooth-time = Yaw reset vertraging (0s schakelt afvlakking uit)
settings-general-tracker_mechanics-drift_compensation = Drift compensatie
# This cares about multilines
settings-general-tracker_mechanics-drift_compensation-description =
    Compenseert voor IMU yaw drift door de toevoeging van een omgekeerde rotatie.
    Veranderd de sterkte van de compensatie en hoeveel resets worden gebruikt.
settings-general-tracker_mechanics-drift_compensation-enabled-label = Drift compensate
settings-general-tracker_mechanics-drift_compensation-prediction = Voorspelling van driftcompensatie
# This cares about multilines
settings-general-tracker_mechanics-drift_compensation-prediction-description =
    Voorspelt compensatie van gierdrift buiten het eerder gemeten bereik.
    Schakel dit in als jouw trackers continu om de gier-as draaien.
settings-general-tracker_mechanics-drift_compensation-prediction-label = Voorspelling van driftcompensatie
settings-general-tracker_mechanics-drift_compensation_warning =
    <b>Waarschuwing:</b> Gebruik alleen driftcompensatie als je heel vaak
    moet resetten (elke ~5-10 minuten).
    
    IMU's die vaak worden gereset, zijn onder ander:
    Joy-Cons, owoTrack en MPU's (zonder recente firmware).
settings-general-tracker_mechanics-drift_compensation_warning-cancel = Annuleren
settings-general-tracker_mechanics-drift_compensation_warning-done = Ik begrijp het
settings-general-tracker_mechanics-drift_compensation-amount-label = Compensatiesterkte
settings-general-tracker_mechanics-drift_compensation-max_resets-label = Gebruik de laatste x resets
settings-general-tracker_mechanics-save_mounting_reset = Sla de automatische montage reset kalibratie op
settings-general-tracker_mechanics-save_mounting_reset-description = Slaat de automatische montage reset kalibraties van de trackers op tussen herstarts. Nuttig als je een pak draagt waarbij trackers niet bewegen tussen sessies. <b>Niet aanbevolen voor normale gebruikers!</b>
settings-general-tracker_mechanics-save_mounting_reset-enabled-label = Montage configuratie opslaan
settings-general-tracker_mechanics-use_mag_on_all_trackers = Gebruik de magnetometer op alle IMU-trackers die dit ondersteunen
settings-general-tracker_mechanics-use_mag_on_all_trackers-description =
    Gebruikt magnetometer op alle trackers die er een compatibele firmware voor hebben, waardoor drift in stabiele magnetische omgevingen wordt verminderd.
    Je kan dit per individuele tracker uit zetten in de instellingen van de tracker. <b>Sluit geen van de trackers af terwijl u dit in- en uitschakelt!</b>
settings-general-tracker_mechanics-use_mag_on_all_trackers-label = Gebruik magnetometer op de trackers
settings-stay_aligned = Blijf in lijn
settings-stay_aligned-description = ijf in lijn vermindert drift door je trackers geleidelijk aan te passen zodat ze overeenkomen met je ontspannen houdingen.
settings-stay_aligned-setup-label = Blijf in lijn instellen
settings-stay_aligned-setup-description = Je moet "Blijf in lijn instellen" voltooien om Blijf in lijn te activeren.
settings-stay_aligned-warnings-drift_compensation = ⚠ Schakel Drift Compensation uit! Drift Compensation conflicteert met Blijf in lijn.
settings-stay_aligned-enabled-label = Trackers aanpassen
settings-stay_aligned-hide_yaw_correction-label = Aanpassing verbergen (om te vergelijken zonder Blijf in lijn)
settings-stay_aligned-general-label = Algemeen
settings-stay_aligned-relaxed_poses-label = Ontspannen houdingen
settings-stay_aligned-relaxed_poses-description = Blijf in lijn gebruikt je ontspannen houdingen om je trackers in lijn te houden. Gebruik "Stel Blijf in lijn in" om deze houdingen bij te werken.
settings-stay_aligned-relaxed_poses-standing = Pas trackers aan terwijl je staat
settings-stay_aligned-relaxed_poses-sitting = Pas trackers aan terwijl je op een stoel zit
settings-stay_aligned-relaxed_poses-flat = Pas trackers aan terwijl je op de grond zit of op je rug ligt.
settings-stay_aligned-relaxed_poses-save_pose = Sla houding op
settings-stay_aligned-relaxed_poses-reset_pose = Reset houding
settings-stay_aligned-relaxed_poses-close = Sluiten
settings-stay_aligned-debug-label = Foutopsporing
settings-stay_aligned-debug-description = Voeg je instellingen toe wanneer je problemen met Blijf in lijn rapporteert.
settings-stay_aligned-debug-copy-label = Instellingen naar klembord kopiëren

## FK/Tracking settings

settings-general-fk_settings = Tracking instellingen
# Floor clip:
# why the name - came from the idea of noclip in video games, but is the opposite where clipping to the floor is a desired feature
# definition - Prevents the foot trackers from going lower than they where when a reset was performed
settings-general-fk_settings-leg_tweak-floor_clip = Floor-clip
# Skating correction:
# why the name - without this enabled the feet will often slide across the ground as if your skating across the ground,
# since this largely prevents this it corrects for it hence skating correction (note this may be renamed to sliding correction)
# definition - Guesses when each foot is in contact with the ground and uses that information to improve tracking
settings-general-fk_settings-leg_tweak-skating_correction = Skating-correctie
settings-general-fk_settings-leg_tweak-toe_snap = Teen snap
settings-general-fk_settings-leg_tweak-foot_plant = Voetplant
settings-general-fk_settings-leg_tweak-skating_correction-amount = Skating-correctie sterkte
settings-general-fk_settings-leg_tweak-skating_correction-description = Schaatscorrectie corrigeert voor schaatsen, maar kan de nauwkeurigheid van bepaalde bewegingspatronen verminderen. Zorg ervoor dat je bij het inschakelen een volledige reset uitvoert en opnieuw kalibreert in het spel.
settings-general-fk_settings-leg_tweak-floor_clip-description =
    Floor-clip kan het doorknippen van de vloer verminderen of zelfs elimineren.
    Zorg ervoor dat u bij het inschakelen een volledige reset uitvoert en opnieuw kalibreert in het spel.
settings-general-fk_settings-leg_tweak-toe_snap-description = Toe-snap probeert de rotatie van uw voeten te raden als voet-trackers niet worden gebruikt.
settings-general-fk_settings-leg_tweak-foot_plant-description = Foot-plant roteert je voeten zodat ze evenwijdig aan de grond zijn wanneer ze in contact zijn.
settings-general-fk_settings-leg_fk = Been tracking
settings-general-fk_settings-leg_fk-reset_mounting_feet-description = Schakel Montage Reset voor de voeten in door op je tenen te staan.
settings-general-fk_settings-leg_fk-reset_mounting_feet = Voeten montage reset.
settings-general-fk_settings-enforce_joint_constraints = Bewegingslimieten van het skelet
settings-general-fk_settings-enforce_joint_constraints-enforce_constraints = Beperkingen toepassen
settings-general-fk_settings-enforce_joint_constraints-enforce_constraints-description = Voorkomt dat gewrichten over hun limiet draaien
settings-general-fk_settings-enforce_joint_constraints-correct_constraints = Corrigeren met beperkingen
settings-general-fk_settings-enforce_joint_constraints-correct_constraints-description = Corrigeer gewrichtsrotaties wanneer ze hun limiet overschrijden
settings-general-fk_settings-arm_fk = Arm tracking
settings-general-fk_settings-arm_fk-description = Verander de manier waarop de armen worden getrackt.
settings-general-fk_settings-arm_fk-force_arms = Dwing armen vanuit HMD
settings-general-fk_settings-reset_settings = Instellingen resetten
settings-general-fk_settings-reset_settings-reset_hmd_pitch-description = Reset de pitch (verticale rotatie) van de HMD na een volledige reset. Dit is handig als je de HMD op je voorhoofd draagt voor VTubing of mocap. Niet inschakelen voor VR.
settings-general-fk_settings-reset_settings-reset_hmd_pitch = HMD pitch resetten
settings-general-fk_settings-arm_fk-reset_mode-description = Pas de verwachte armhouding aan voor het resetten van de montage.
settings-general-fk_settings-arm_fk-back = Achterzijde
settings-general-fk_settings-arm_fk-back-description = De standaardmodus, waarbij de bovenarmen  naar achteren gaan en de onderarmen naar voren.
settings-general-fk_settings-arm_fk-tpose_up = T-pose (omhoog)
settings-general-fk_settings-arm_fk-tpose_up-description = Verwacht je armen langs je zeiden te hangen tijdens een volledige reset, en 90 graden omhoog langs je zeiden tijdens een Montage Reset.
settings-general-fk_settings-arm_fk-tpose_down = T-pose (omlaag)
settings-general-fk_settings-arm_fk-tpose_down-description = Verwacht dat je armen 90 graden naar de zijkanten zijn tijdens een Volledige reset, en aan de zijkanten naar beneden tijdens een montage reset.
settings-general-fk_settings-arm_fk-forward = Voorwaards
settings-general-fk_settings-arm_fk-forward-description = Verwacht dat je armen 90 graden naar voren staan. Handig voor VTubing.
settings-general-fk_settings-skeleton_settings-toggles = Skelet schakelaars
settings-general-fk_settings-skeleton_settings-description = Schakel skeleton instellingen in of uit. Het is aanbevolen om deze aan te laten.
settings-general-fk_settings-skeleton_settings-extended_spine_model = Uitgebreid ruggengraat model
settings-general-fk_settings-skeleton_settings-extended_pelvis_model = Uitgebreid bekken model
settings-general-fk_settings-skeleton_settings-extended_knees_model = Uitgebreid knie model
settings-general-fk_settings-skeleton_settings-ratios = skelet verhoudingen
settings-general-fk_settings-skeleton_settings-ratios-description = Pas de waardes van de skelet instellingen aan. Het kan zijn dat je hierna je lichaams proporties moet aanpassen.
settings-general-fk_settings-skeleton_settings-impute_waist_from_chest_hip = Bereken taille van borst naar heup
settings-general-fk_settings-skeleton_settings-impute_waist_from_chest_legs = Bereken taille van borst naar benen
settings-general-fk_settings-skeleton_settings-impute_hip_from_chest_legs = Bereken heup van borst naar benen
settings-general-fk_settings-skeleton_settings-impute_hip_from_waist_legs = Bereken heup van taille naar benen
settings-general-fk_settings-skeleton_settings-interp_hip_legs = Bereken het gemiddelde van de 'yaw en roll van de heup met die van de benen'
settings-general-fk_settings-skeleton_settings-interp_knee_tracker_ankle = Bereken het gemiddelde van de 'yaw en roll van de knie trackers met die van de enkels'
settings-general-fk_settings-skeleton_settings-interp_knee_ankle = Bereken het gemiddelde van de 'yaw en roll van de knie trackers met die van de enkels'
settings-general-fk_settings-self_localization-title = Mocap modus
settings-general-fk_settings-self_localization-description = Mocap modus staat het skelet model toe om zijn eigen positie te bepalen zonder het gebruik van een headset of andere trackers. Dit vergt wel het gebruik van voet en hoofd trackers, dit is momenteel nog expirimenteel.

## Gesture control settings (tracker tapping)

settings-general-gesture_control = Gesture control
settings-general-gesture_control-subtitle = Op tik gebaseerde resets
settings-general-gesture_control-description = Maakt het mogelijk om resets te activeren door op een tracker te tikken. De tracker het hoogst op je bovenlichaam wordt gebruikt voor Quick Reset, de tracker het hoogst op je linkerbeen voor Reset en de tracker het hoogst op je rechterbeen voor Montage Reset. Het moet worden vermeld dat tikken binnen 0,6 seconden moeten gebeuren om geregistreerd te worden.
# This is a unit: 3 taps, 2 taps, 1 tap
# $amount (Number) - Amount of taps (touches to the tracker's case)
settings-general-gesture_control-taps =
    { $amount ->
        [one] 1 tik
       *[other] { $amount } tikken
    }
# This is a unit: 3 trackers, 2 trackers, 1 tracker
# $amount (Number) - Amount of trackers
settings-general-gesture_control-trackers =
    { $amount ->
        [one] één
       *[other] anders
    }
settings-general-gesture_control-yawResetEnabled = Activeer tikken voor horizontale reset
settings-general-gesture_control-yawResetDelay = Vertraging horizontale reset
settings-general-gesture_control-yawResetTaps = Hoeveelheid tikken voor horizontale reset
settings-general-gesture_control-fullResetEnabled = Activeer tikken voor volledige reset
settings-general-gesture_control-fullResetDelay = Vertraging volledige reset
settings-general-gesture_control-fullResetTaps = Hoeveelheid tikken voor volledige reset
settings-general-gesture_control-mountingResetEnabled = Activeer tikken voor montage-kalibratie
settings-general-gesture_control-mountingResetDelay = Vertraging montage-kalibratie
settings-general-gesture_control-mountingResetTaps = Hoeveelheid tikken voor montage-kalibratie
# The number of trackers that can have higher acceleration before a tap is rejected
settings-general-gesture_control-numberTrackersOverThreshold = Trackers over drempelwaarde
settings-general-gesture_control-numberTrackersOverThreshold-description = Verhoog deze waarde als de tik detectie niet werkt. Zet deze waarde niet te hoog om tik detectie te laten werken, dit kan vals positieve resultaten creëren.

## Appearance settings

settings-interface-appearance = Uiterlijk
settings-general-interface-dev_mode = Ontwikkelaarsmodus
settings-general-interface-dev_mode-description = Deze modus kan nuttig zijn als je diepgaande gegevens nodig hebt of op een geavanceerd niveau wilt communiceren met aangesloten trackers.
settings-general-interface-dev_mode-label = Ontwikkelaarsmodus
settings-general-interface-theme = Themakleur
settings-general-interface-show-navbar-onboarding = Toon "{ navbar-onboarding }" op de navigatiebalk
settings-general-interface-show-navbar-onboarding-description = Dit verandert of de knop "{ navbar-onboarding }" wordt weergegeven op de navigatiebalk.
settings-general-interface-show-navbar-onboarding-label = Toon "{ navbar-onboarding }"
settings-general-interface-lang = Selecteer taal
settings-general-interface-lang-description = Verander de standaardtaal die je wilt gebruiken.
settings-general-interface-lang-placeholder = Selecteer de te gebruiken taal
# Keep the font name untranslated
settings-interface-appearance-font = GUI lettertype
settings-interface-appearance-font-description = Dit past het lettertype aan welke gebruikt wordt binnen het interface
settings-interface-appearance-font-placeholder = Standaard lettertype
settings-interface-appearance-font-os_font = Besturingssysteem lettertype
settings-interface-appearance-font-slime_font = Standaard lettertype
settings-interface-appearance-font_size = Standaard lettertype grote
settings-interface-appearance-font_size-description = Dit past het lettertype grote aan voor het gehele interfeace, behalve voor deze instellingen pagina.
settings-interface-appearance-decorations = Gebruik het systeem native decoraties
settings-interface-appearance-decorations-description = Dit zal de bovenste balk van de interface niet weergeven en zal in plaats daarvan die van het besturingssysteem gebruiken.
settings-interface-appearance-decorations-label = Gebruik de native decoraties

## Notification settings

settings-interface-notifications = Notificaties
settings-general-interface-serial_detection = Detectie van seriële apparaten
settings-general-interface-serial_detection-description = Met deze optie verschijnt er elke keer dat je een nieuw serieel apparaat aansluit dat mogelijk een tracker is, een pop-up. Dit helpt bij het verbeteren van het instelproces van een tracker.
settings-general-interface-serial_detection-label = Detectie van seriële apparaten
settings-general-interface-feedback_sound = Feedback geluid
settings-general-interface-feedback_sound-description = Speelt een geluid telkens de reset wordt uitgevoerd
settings-general-interface-feedback_sound-label = Feedback geluid
settings-general-interface-feedback_sound-volume = Feedback geluid volume
settings-general-interface-connected_trackers_warning = Waarschuwing voor verbonden trackers
settings-general-interface-connected_trackers_warning-description = Deze optie toont een pop-up bericht telkens wanneer je SlimeVR probeert af te sluiten terwijl er nog trackers verbonden zijn. Dit bericht herinnert je eraan om je trackers uit te schakelen wanneer je klaar bent om de batterijduur te sparen.
settings-general-interface-connected_trackers_warning-label = Waarschuwing voor verbonden trackers bij het afsluiten

## Behavior settings

settings-interface-behavior = Gedrag
settings-general-interface-dev_mode = Ontwikkelaarsmodus
settings-general-interface-dev_mode-description = Deze modus kan nuttig zijn als je diepgaande gegevens nodig hebt of op een geavanceerd niveau wilt communiceren met aangesloten trackers.
settings-general-interface-dev_mode-label = Ontwikkelaarsmodus
settings-general-interface-use_tray = Minimaliseren naar systeem vak
settings-general-interface-use_tray-description = Hiermee kun je het venster sluiten zonder de SlimeVR server te beëindigen, zodat je deze op de achtergrond kunt blijven gebruiken zonder dat de GUI in de weg zit.
settings-general-interface-use_tray-label = Minimaliseren naar systeem vak
settings-general-interface-discord_presence = Activiteit delen op Discord
settings-general-interface-discord_presence-description = Deelt op Discord dat je de SlimeVR server gebruikt, tezamen met het aantal IMU-Trackers.
settings-general-interface-discord_presence-label = Activiteit delen op Discord
settings-general-interface-discord_presence-message =
    { $amount ->
        [0] Aan het slimen
        [one] Gebruikt 1 tracker
       *[other] Gebruikt { $amount } trackers
    }
settings-interface-behavior-error_tracking = Foutverzameling via Sentry.io
settings-interface-behavior-error_tracking-description_v2 =
    <h1>Geeft u toestemming voor het verzamelen van geanonimiseerde foutgegevens?</h1>
    
    <b>We verzamelen geen persoonlijke informatie</b> zoals uw IP-adres of draadloze inloggegevens. SlimeVR hecht veel waarde aan uw privacy!
    
    Om de beste gebruikerservaring te bieden, verzamelen we geanonimiseerde foutrapporten, prestatiestatistieken en informatie over het besturingssysteem. Dit helpt ons bij het detecteren van fouten en problemen met SlimeVR. Deze statistieken worden verzameld via Sentry.io.
settings-interface-behavior-error_tracking-label = Stuur fouten naar de ontwikkelaars
settings-interface-behavior-bvh_directory = Map om BVH-opnames op te slaan
settings-interface-behavior-bvh_directory-description = Kies een map om je BVH-opnames op te slaan, zodat je niet elke keer hoeft te kiezen waar je ze opslaat.
settings-interface-behavior-bvh_directory-label = Map voor BVH-opnames

## Serial settings

settings-serial = Seriele console
# This cares about multilines
settings-serial-description =
    Dit is een live informatiefeed voor seriële communicatie.
    Kan handig zijn als je wilt weten dat de firmware werkt.
settings-serial-connection_lost = Verbinding met seriële poort verloren, opnieuw verbinden...
settings-serial-reboot = Opnieuw opstarten
settings-serial-factory_reset = Fabrieksinstellingen herstellen
# This cares about multilines
# <b>text</b> means that the text should be bold
settings-serial-factory_reset-warning =
    <b>Waarschuwing:</b> Hiermee wordt de tracker teruggezet naar de fabrieksinstellingen.
    Wat betekent dat Wi-Fi en kalibratie-instellingen <b>allemaal verloren gaan!</b>
settings-serial-factory_reset-warning-ok = Ik weet wat ik doe
settings-serial-factory_reset-warning-cancel = Annuleren
settings-serial-get_infos = Informatie ophalen
settings-serial-serial_select = Selecteer een seriële poort
settings-serial-auto_dropdown_item = Automatisch
settings-serial-get_wifi_scan = WiFi-scan uitvoeren
settings-serial-file_type = Gewone tekst
settings-serial-save_logs = Opslaan in bestand

## OSC router settings

settings-osc-router = OSC-router
# This cares about multilines
settings-osc-router-description =
    Stuur OSC-berichten door vanuit een ander programma.
    Nuttig om bijvoorbeeld een ander OSC-programma te gebruiken met VRChat.
settings-osc-router-enable = Inschakelen
settings-osc-router-enable-description = Schakel het doorsturen van berichten in of uit.
settings-osc-router-enable-label = Inschakelen
settings-osc-router-network = Netwerkpoorten
# This cares about multilines
settings-osc-router-network-description =
    Stel de poorten in voor het verzenden en ontvangen van gegevens.
    Dit kunnen dezelfde poorten zijn als andere poorten die worden gebruikt in de SlimeVR-server.
settings-osc-router-network-port_in =
    .label = Poort in
    .placeholder = Poort in (standaard: 9002)
settings-osc-router-network-port_out =
    .label = Poort uit
    .placeholder = Poort uit (standaard: 9000)
settings-osc-router-network-address = Netwerkadres
settings-osc-router-network-address-description = Stel het adres in waarnaar gegevens moeten worden verzonden.
settings-osc-router-network-address-placeholder = IPV4-adres

## OSC VRChat settings

settings-osc-vrchat = VRChat OSC Trackers
# This cares about multilines
settings-osc-vrchat-description-v1 = Wijzig instellingen die specifiek zijn voor de OSC Trackers-standaard die wordt gebruikt voor het verzenden van trackinggegevens naar applicaties zonder SteamVR (bijv. Quest standalone). Zorg ervoor dat OSC is ingeschakeld in VRChat via het Actiemenu onder OSC > Ingeschakeld. Om het ontvangen van HMD- en controllergegevens van VRChat mogelijk te maken, ga in je hoofdmenu naar Instellingen onder Tracking & IK > Allow Sending Head and Wrist VR Tracking OSC Data.
settings-osc-vrchat-enable = Inschakelen
settings-osc-vrchat-enable-description = Schakel het verzenden en ontvangen van gegevens in en uit.
settings-osc-vrchat-enable-label = Inschakelen
settings-osc-vrchat-oscqueryEnabled = OSCQuery inschakelen
settings-osc-vrchat-oscqueryEnabled-description =
    OSCQuery detecteert automatisch actieve instanties van VRChat en stuurt hen data.
    Het kan zichzelf ook aan hen bekendmaken om HMD- en controllerdata te ontvangen.
    Om HMD- en controllerdata van VRChat te kunnen ontvangen, ga je in het hoofdmenu naar Instellingen,
    onder "Tracking & IK", en schakel je "Allow Sending Head and Wrist VR Tracking OSC Data" in.
settings-osc-vrchat-oscqueryEnabled-label = Schakel OSCQuery in
settings-osc-vrchat-network = Netwerkpoorten
settings-osc-vrchat-network-description-v1 = Stel de poorten in voor het ontvangen en verzenden van tracking data. Kan op standaardinstellingen blijven voor VRChat.
settings-osc-vrchat-network-port_in =
    .label = Poort In
    .placeholder = Poort in (standaard: 9001)
settings-osc-vrchat-network-port_out =
    .label = Poort Out
    .placeholder = Poort uit (standaard: 9000)
settings-osc-vrchat-network-address = Netwerkadres
settings-osc-vrchat-network-address-description-v1 = Kies naar welk adres u gegevens wilt verzenden. Kan op standaardinstellingen blijven voor VRChat.
settings-osc-vrchat-network-address-placeholder = VRChat IP-adres
settings-osc-vrchat-network-trackers = Trackers
settings-osc-vrchat-network-trackers-description = Schakel het verzenden van specifieke trackers via OSC in en uit.
settings-osc-vrchat-network-trackers-chest = Borst
settings-osc-vrchat-network-trackers-hip = Heup
settings-osc-vrchat-network-trackers-knees = Knieën
settings-osc-vrchat-network-trackers-feet = Voeten
settings-osc-vrchat-network-trackers-elbows = Ellebogen

## VMC OSC settings

settings-osc-vmc = Virtuele motion capture
# This cares about multilines
settings-osc-vmc-description =
    Verander instellingen specifiek voor het VMC (Virtual Motion Capture) protocol
     botgegevens van SlimeVR te verzenden en botgegevens van andere apps te ontvangen.
settings-osc-vmc-enable = Inschakelen
settings-osc-vmc-enable-description = Schakel het verzenden en ontvangen van gegevens in en uit.
settings-osc-vmc-enable-label = Inschakelen
settings-osc-vmc-network = Netwerkpoorten
settings-osc-vmc-network-description = Stel de poorten in voor het zenden en ontvangen van VMC-gegevens.
settings-osc-vmc-network-port_in =
    .label = Poort In
    .placeholder = Poort in (standaard: 39540)
settings-osc-vmc-network-port_out =
    .label = Poort uit
    .placeholder = Poort uit (standaard: 39539)
settings-osc-vmc-network-address = Netwerkadres
settings-osc-vmc-network-address-description = Stel het adres in waarnaar gegevens moeten worden verzonden via VMC.
settings-osc-vmc-network-address-placeholder = IPV4-adres
settings-osc-vmc-vrm = VRM Model
settings-osc-vmc-vrm-description = Laad een VRM-model om hoofdverankering mogelijk te maken en zorg voor een hogere compatibiliteit met andere applicaties.
settings-osc-vmc-vrm-untitled_model = Naamloos model
settings-osc-vmc-vrm-file_select = Sleep een modelbestand naar hier om ze te gebruiken of <u>blader</u>.
settings-osc-vmc-anchor_hip = Heupverankering
settings-osc-vmc-anchor_hip-description = Veranker de tracking aan de heupen, handig voor zittende VTubing. Als u deze uitschakelt, laadt u een VRM-model.
settings-osc-vmc-anchor_hip-label = Heupverankering
settings-osc-vmc-mirror_tracking = Gespiegelde tracking
settings-osc-vmc-mirror_tracking-description = De tracking horizontaal spiegelen.
settings-osc-vmc-mirror_tracking-label = Gespiegelde tracking

## Advanced settings

settings-utils-advanced = Geavanceerd
settings-utils-advanced-reset-gui = GUI-instellingen resetten
settings-utils-advanced-reset-gui-description = Herstel de standaardwaarden voor instellingen van de interface.
settings-utils-advanced-reset-gui-label = GUI resetten
settings-utils-advanced-reset-server = Tracker instellingen resetten
settings-utils-advanced-reset-server-description = Herstel de standaardwaarden voor instellingen van de tracking.
settings-utils-advanced-reset-server-label = Tracking resetten
settings-utils-advanced-reset-all = Alle instellingen resetten
settings-utils-advanced-reset-all-description = Herstel de standaardwaarden voor instellingen van beide de interface en de tracking.
settings-utils-advanced-reset-all-label = Alles resetten
settings-utils-advanced-reset_warning =
    { $type ->
        [gui]
            <b>Waarschuwing</b>Hiermee worden al uw GUI instellingen teruggezet naar de standaardinstellingen. 
            Weet u zeker dat u dit wilt doen?
        [server]
            <b>Waarschuwing</b>Hiermee worden al uw tracking instellingen teruggezet naar de standaardinstellingen.
            Weet u zeker dat u dit wilt doen?
       *[all]
            <b>Waarschuwing:</b> Hiermee worden al uw instellingen teruggezet naar de standaardinstellingen.
            Weet u zeker dat u dit wilt doen?
    }
settings-utils-advanced-reset_warning-reset = Instellingen resetten
settings-utils-advanced-reset_warning-cancel = Annuleren
settings-utils-advanced-open_data-v1 = Configuratiemap
settings-utils-advanced-open_data-description-v1 = Open de configuratiemap van SlimeVR in de bestandsverkenner, met configuratiebestanden.
settings-utils-advanced-open_data-label = Map openen
settings-utils-advanced-open_logs = logboeken
settings-utils-advanced-open_logs-description = Open de logmap van SlimeVR in de bestandsverkenner, met de logboeken van de app
settings-utils-advanced-open_logs-label = Map openen

## Setup/onboarding menu

onboarding-skip = Setupgids overslaan
onboarding-continue = Doorgaan
onboarding-wip = WIP
onboarding-previous_step = Vorige stap
onboarding-setup_warning =
    <b>Waarschuwing:</b> De initiële setup is nodig voor een goede tracking ervaring,
    het is aangeraden deze te volgen indien dit de eerste keer is dat je SlimeVR gebruikt.
onboarding-setup_warning-skip = Setupgids overslaan
onboarding-setup_warning-cancel = Doorgaan met setupgids

## Wi-Fi setup

onboarding-wifi_creds-back = Ga terug naar de introductie
onboarding-wifi_creds = Voer de WiFi-inloggegevens in
# This cares about multilines
onboarding-wifi_creds-description =
    Deze gegevens worden gebruikt om de trackers draadloos te verbinden met de server.
    Gelieve de gegevens te gebruiken van het netwerk waarmee je momenteel bent verbonden.
onboarding-wifi_creds-skip = WiFi-instellingen overslaan
onboarding-wifi_creds-submit = Verzenden!
onboarding-wifi_creds-ssid =
    .label = WiFi naam
    .placeholder = Vul WiFi naam in
onboarding-wifi_creds-ssid-required = Wi-Fi-naam is vereist
onboarding-wifi_creds-password =
    .label = Paswoord
    .placeholder = Vul paswoord in

## Mounting setup

onboarding-reset_tutorial-back = Ga terug naar de montage-kalibratie
onboarding-reset_tutorial = Reset tutorial
onboarding-reset_tutorial-explanation = Terwijl je jouw trackers gebruikt, kunnen ze uit de lijn raken vanwege IMU-yaw-drift, of omdat je ze fysiek hebt verplaatst. Je hebt verschillende manieren om dit op te lossen.
onboarding-reset_tutorial-skip = Stap overslaan
# Cares about multiline
onboarding-reset_tutorial-0 =
    Tik { $taps } keer op de gemarkeerde tracker om de yaw-reset te activeren.
    
    Hierdoor staan de trackers in dezelfde richting als je HMD.
# Cares about multiline
onboarding-reset_tutorial-1 =
    Tik { $taps } keer op de gemarkeerde tracker om een volledige reset uit te voeren.
    
    Hiervoor moet je staan (i-pose). Er is een vertraging van 3 seconden (instelbaar) voordat het daadwerkelijk gebeurt.
    Hiermee wordt de positie en rotatie van al je trackers volledig gereset. Dit zou de meeste problemen moeten oplossen.
# Cares about multiline
onboarding-reset_tutorial-2 =
    Tik { $taps } keer op de gemarkeerde tracker om de montage opnieuw in te stellen.
    
    Montage-reset helpt bij hoe de trackers daadwerkelijk op je worden geplaatst, dus als je ze per ongeluk hebt verplaatst en de oriëntatie ervan voor een groot deel hebt veranderd, zal dit helpen.
    
    Je moet in een houding staan alsof je aan het skiën bent, zoals wordt weergegeven in de Automatische montage wizard, je hebt een vertraging van 3 seconden (instelbaar) voordat deze wordt geactiveerd.

## Setup start

onboarding-home = Welkom bij SlimeVR
onboarding-home-start = Laten we beginnen!

## Enter VR part of setup

onboarding-enter_vr-back = Ga terug naar de sectie voor toewijzing van trackers
onboarding-enter_vr-title = Tijd om VR in te gaan!
onboarding-enter_vr-description = Doe al je trackers aan en ga dan in VR!
onboarding-enter_vr-ready = Gereed!

## Setup done

onboarding-done-title = Je bent klaar!
onboarding-done-description = Geniet van je full-body ervaring
onboarding-done-close = Sluit de gids

## Tracker connection setup

onboarding-connect_tracker-back = Ga terug naar de instellingen voor WiFi-configuratie
onboarding-connect_tracker-title = Trackers verbinden
onboarding-connect_tracker-description-p0-v1 = Op naar het leukste gedeelte, trackers verbinden!
onboarding-connect_tracker-description-p1-v1 = Sluit elke tracker één voor één aan via een USB poort.
onboarding-connect_tracker-issue-serial = Ik heb problemen met verbinden!
onboarding-connect_tracker-usb = USB Tracker
onboarding-connect_tracker-connection_status-none = Op zoek naar trackers
onboarding-connect_tracker-connection_status-serial_init = Verbinding maken met een serieel apparaat
onboarding-connect_tracker-connection_status-obtaining_mac_address = Het mac-adres van de tracker verkrijgen
onboarding-connect_tracker-connection_status-provisioning = Wifi-inloggegevens verzenden
onboarding-connect_tracker-connection_status-connecting = Wifi-inloggegevens verzenden
onboarding-connect_tracker-connection_status-looking_for_server = Op zoek naar server
onboarding-connect_tracker-connection_status-connection_error = Kan geen verbinding maken met Wi-Fi
onboarding-connect_tracker-connection_status-could_not_find_server = Kan de server niet vinden
onboarding-connect_tracker-connection_status-done = Verbonden met de server
onboarding-connect_tracker-connection_status-no_serial_log = Kon geen logbestanden van de tracker ophalen.
onboarding-connect_tracker-connection_status-no_serial_device_found = Kon geen tracker via USB vinden.
onboarding-connect_serial-error-modal-no_serial_log = Staat de tracker aan?
onboarding-connect_serial-error-modal-no_serial_log-desc = Zorg dat de tracker aan staat en verbonden is met je computer.
onboarding-connect_serial-error-modal-no_serial_device_found = Geen trackers gedetecteerd
onboarding-connect_serial-error-modal-no_serial_device_found-desc =
    Sluit een tracker met de meegeleverde USB-kabel aan op je computer en zet de tracker aan.
    Als dit niet werkt:
      -probeer een andere USB-kabel
      -probeer een andere USB-poort
      -probeer de SlimeVR-server opnieuw te installeren en selecteer "USB Drivers" in het onderdeelkeuze-menu
# $amount (Number) - Amount of trackers connected (this is a number, but you can use CLDR plural rules for your language)
# More info on https://www.unicode.org/cldr/cldr-aux/charts/22/supplemental/language_plural_rules.html
# English in this case only has 2 plural rules, which are "one" and "other",
# we use 0 in an explicit way because there is no plural rule in english for 0, so we directly say
# if $amount is 0 then we say "No trackers connected"
onboarding-connect_tracker-connected_trackers =
    { $amount ->
        [0] Geen trackers
        [one] 1 tracker
       *[other] { $amount } trackers
    } verbonden
onboarding-connect_tracker-next = Ik heb al mijn trackers verbonden

## Tracker calibration tutorial

onboarding-calibration_tutorial = Handleiding voor IMU-kalibratie
onboarding-calibration_tutorial-subtitle = Helpt met het verminderen van het driften van de trackers!
onboarding-calibration_tutorial-description-v1 = Zet je trackers aan en leg ze even op een stabiele ondergrond om te kalibreren. Kalibratie kan op elk moment na het inschakelen van de trackers worden uitgevoerd — deze pagina biedt alleen een handleiding. Klik op de knop "{ onboarding-calibration_tutorial-calibrate }" om te beginnen en <b>beweeg je trackers daarna niet meer!</b>
onboarding-calibration_tutorial-calibrate = Al mijn trackers liggen neer
onboarding-calibration_tutorial-status-waiting = Ik wacht op jou
onboarding-calibration_tutorial-status-calibrating = Kalibreren
onboarding-calibration_tutorial-status-success = Aardig!
onboarding-calibration_tutorial-status-error = De tracker werd verplaatst
onboarding-calibration_tutorial-skip = Tutorial overslaan

## Tracker assignment tutorial

onboarding-assignment_tutorial = Hoe een Slime Tracker voor te bereiden voordat u deze aantrekt
onboarding-assignment_tutorial-first_step = 1. Plaats een lichaamsdeelsticker (als je die hebt) op de tracker naar keuze
# This text has a character limit of around 11 characters, so please keep it short
onboarding-assignment_tutorial-sticker = Sticker
onboarding-assignment_tutorial-second_step-v2 = Bevestig de strap aan de tracker met de klittenbandzijde in dezelfde richting als de voorzijde van de tracker:
onboarding-assignment_tutorial-second_step-continuation-v2 = De klittenbandzijde van de extensie moet naar boven gericht zijn, zoals in de foto afgebeeld:
onboarding-assignment_tutorial-done = Ik heb stickers en riemen geplaatst!

## Tracker assignment setup

onboarding-assign_trackers-back = Ga terug naar de instellingen voor WiFi-configuratie
onboarding-assign_trackers-title = Trackers toewijzen
onboarding-assign_trackers-description = Laten we de bevesteging van je trackers bepalen. Klik op de lichaamslocatie waar je een tracker wilt toewijzen.
# Look at translation of onboarding-connect_tracker-connected_trackers on how to use plurals
# $assigned (Number) - Trackers that have been assigned a body part
# $trackers (Number) - Trackers connected to the server
onboarding-assign_trackers-assigned =
    { $assigned } van { $trackers ->
        [one] 1 tracker
       *[other] { $trackers } trackers
    } toegewezen
onboarding-assign_trackers-advanced = Geavanceerde toewijzingslocaties weergeven
onboarding-assign_trackers-next = Ik heb alle trackers toegewezen
onboarding-assign_trackers-mirror_view = Gespiegelde weergave
onboarding-assign_trackers-option-amount =
    { $trackersCount ->
        [one] x{ $trackersCount }
       *[other] x{ $trackersCount }
    }
onboarding-assign_trackers-option-label =
    { $mode ->
        [lower-body] Lower-Body Set
        [core] Core Set
        [enhanced-core] Enhanced Core Set
        [full-body] Full-Body Set
       *[all] Alle trackers
    }
onboarding-assign_trackers-option-description =
    { $mode ->
        [lower-body] Minimaal aantal voor VR full-body tracking
        [core] + betere torso tracking
        [enhanced-core] + voeten rotatie
        [full-body] + elleboog tracking
       *[all] Alle beschikbare tracker locaties
    }

## Tracker assignment warnings

# Note for devs, number is used for representing boolean states per bit.
# $unassigned (Number) - Bits are based on BodyAssignment.ASSIGNMENT_RULES order
onboarding-assign_trackers-warning-LEFT_FOOT =
    { $unassigned ->
        [0] De linkervoet is toegewezen, maar de linkerenkel, linkerdij en de borst, heup of taille moeten ook worden toegewezen!
        [1] De linkervoet is toegewezen, maar het linkerdij en de borst, heup of taille moeten ook worden toegewezen!
        [2] De linkervoet is toegewezen, maar de linkerenkel en de borst, heup of taille moeten ook worden toegewezen!
        [3] De linkervoet is toegewezen, maar de borst, heup of taille moeten ook worden toegewezen!
        [4] De linkervoet is toegewezen, maar de linkerenkel en linkerdij moeten ook worden toegewezen!
        [5] De linkervoet is toegewezen, maar het linkerdij moet ook worden toegewezen!
        [6] De linkervoet is toegewezen, maar de linkerenkel moet ook worden toegewezen!
       *[other] De linkervoet is toegewezen, maar het onbekend lichaamsdeel moet ook worden toegewezen!
    }
# $unassigned (Number) - Bits are based on BodyAssignment.ASSIGNMENT_RULES order
onboarding-assign_trackers-warning-RIGHT_FOOT =
    { $unassigned ->
        [0] De rechtervoet is toegewezen, maar de rechterenkel, rechterdij en de borst, heup of taille moeten ook worden toegewezen!
        [1] De rechtervoet is toegewezen, maar het rechterdij en de borst, heup of taille moeten ook worden toegewezen!
        [2] De rechtervoet is toegewezen, maar de rechterenkel en de borst, heup of taille moeten ook worden toegewezen!
        [3] De rechtervoet is toegewezen, maar de borst, heup of taille moeten ook worden toegewezen!
        [4] De rechtervoet is toegewezen, maar de rechterenkel en rechterdij moeten ook worden toegewezen!
        [5] De rechtervoet is toegewezen, maar het rechterdij moet ook worden toegewezen!
        [6] De rechtervoet is toegewezen, maar de rechterenkel moet ook worden toegewezen!
       *[other] De rechtervoet is toegewezen, maar het onbekend lichaamsdeel moet ook worden toegewezen!
    }
# $unassigned (Number) - Bits are based on BodyAssignment.ASSIGNMENT_RULES order
onboarding-assign_trackers-warning-LEFT_LOWER_LEG =
    { $unassigned ->
        [0] De linkerenkel is toegewezen, maar de linkerdij en de borst, heup of taille moeten ook worden toegewezen!
        [1] De linkerenkel is toegewezen, maar de borst, heup of taille moeten ook worden toegewezen!
        [2] De linkerenkel is toegewezen, maar de linkerdij moet ook worden toegewezen!
       *[other] De linkerenkel is toegewezen, maar het onbekend lichaamsdeel moet ook worden toegewezen!
    }
# $unassigned (Number) - Bits are based on BodyAssignment.ASSIGNMENT_RULES order
onboarding-assign_trackers-warning-RIGHT_LOWER_LEG =
    { $unassigned ->
        [0] De rechterenkel is toegewezen, maar de rechterdij en de borst, heup of taille moeten ook worden toegewezen!
        [1] De rechterenkel is toegewezen, maar de borst, heup of taille moeten ook worden toegewezen!
        [2] De rechterenkel is toegewezen, maar de rechterdij moet ook worden toegewezen!
       *[other] De rechterenkel is toegewezen, maar het onbekend lichaamsdeel moet ook worden toegewezen!
    }
# $unassigned (Number) - Bits are based on BodyAssignment.ASSIGNMENT_RULES order
onboarding-assign_trackers-warning-LEFT_UPPER_LEG =
    { $unassigned ->
        [0] De linkerdij is toegewezen, maar de borst, heup of taille moeten ook worden toegewezen!
       *[other] De linkerdij is toegewezen, maar het onbekend lichaamsdeel moet ook worden toegewezen!
    }
# $unassigned (Number) - Bits are based on BodyAssignment.ASSIGNMENT_RULES order
onboarding-assign_trackers-warning-RIGHT_UPPER_LEG =
    { $unassigned ->
        [0] De rechterdij is toegewezen, maar de borst, heup of taille moeten ook worden toegewezen!
       *[other] De rechterdij is toegewezen, maar het onbekend lichaamsdeel moet ook worden toegewezen!
    }
# $unassigned (Number) - Bits are based on BodyAssignment.ASSIGNMENT_RULES order
onboarding-assign_trackers-warning-HIP =
    { $unassigned ->
        [0] De heup is toegewezen, maar de borst moet ook worden toegewezen!
       *[other] De heup is toegewezen, maar het onbekend lichaamsdeel moet ook worden toegewezen!
    }
# $unassigned (Number) - Bits are based on BodyAssignment.ASSIGNMENT_RULES order
onboarding-assign_trackers-warning-WAIST =
    { $unassigned ->
        [0] De taille is toegewezen, maar de borst moet ook worden toegewezen!
       *[other] De taille is toegewezen, maar het onbekend lichaamsdeel moet ook worden toegewezen!
    }

## Tracker mounting method choose

onboarding-choose_mounting = Welke montagekalibratiemethode moet worden gebruikt?
# Multiline text
onboarding-choose_mounting-description = De oriëntatie van de montage corrigeert de plaatsing van trackers op je lichaam.
onboarding-choose_mounting-auto_mounting = Automatische bevestiging
# Italicized text
onboarding-choose_mounting-auto_mounting-label-v2 = Aanbevolen
onboarding-choose_mounting-auto_mounting-description = Dit detecteert automatisch de montagerichtingen voor al uw trackers door middel van 2 poses
onboarding-choose_mounting-manual_mounting = Handmatige bevestiging
# Italicized text
onboarding-choose_mounting-manual_mounting-label-v2 = Misschien niet precies genoeg
onboarding-choose_mounting-manual_mounting-description = Hiermee kunt u de montagerichting handmatig kiezen voor elke tracker
# Multiline text
onboarding-choose_mounting-manual_modal-title = Ben je zeker dat je de automatische kalibratie wilt uitvoeren?
onboarding-choose_mounting-manual_modal-description = <b>De handmatige montagekalibratie word aangeraden voor nieuwe gebruikers.</b>, De posities die je moet doen voor de automatische kalibratie kunnen lastig zijn om in één keer goed te krijgen en vereisen mogelijk wat oefening.
onboarding-choose_mounting-manual_modal-confirm = Ik weet zeker wat ik doe
onboarding-choose_mounting-manual_modal-cancel = Annuleren

## Tracker manual mounting setup

onboarding-manual_mounting-back = Ga terug naar de VR sectie
onboarding-manual_mounting = Handmatige montage
onboarding-manual_mounting-description = Klik op elke tracker en selecteer op welke manier ze zijn bevestigd
onboarding-manual_mounting-auto_mounting = Automatische montage
onboarding-manual_mounting-next = Volgende stap

## Tracker automatic mounting setup

onboarding-automatic_mounting-back = Ga terug naar de VR sectie
onboarding-automatic_mounting-title = Montage-kalibratie
onboarding-automatic_mounting-description = Om je trackers te laten werken, moet de rotatie worden ingesteld hoe deze zijn bevestigd op je lichaam.
onboarding-automatic_mounting-manual_mounting = Montage handmatig instellen
onboarding-automatic_mounting-next = Volgende stap
onboarding-automatic_mounting-prev_step = Vorige stap
onboarding-automatic_mounting-done-title = Montagerichtingen gekalibreerd.
onboarding-automatic_mounting-done-description = Je montage-kalibratie is compleet!
onboarding-automatic_mounting-done-restart = Terug naar start
onboarding-automatic_mounting-mounting_reset-title = Montage-reset
onboarding-automatic_mounting-mounting_reset-step-0 = 1. Ga staan in een "skie"-houding met gebogen benen, je bovenlichaam naar voren gekanteld en armen gebogen.
onboarding-automatic_mounting-mounting_reset-step-1 = 2. Druk op de knop "Reset montage" en wacht 3 seconden voordat de montagerichtingen van de trackers opnieuw worden ingesteld.
onboarding-automatic_mounting-preparation-title = Voorbereiding
onboarding-automatic_mounting-preparation-v2-step-0 = 1. Druk op de knop "Volledige reset".
onboarding-automatic_mounting-preparation-v2-step-1 = 2. Ga rechtop staan met je armen langs je zij. Zorg dat je recht vooruit kijkt.
onboarding-automatic_mounting-preparation-v2-step-2 = 3. Houd deze houding aan totdat de timer van 3 seconden is afgelopen.
onboarding-automatic_mounting-put_trackers_on-title = Doe je trackers aan
onboarding-automatic_mounting-put_trackers_on-description = Om montagerichtingen te kalibreren gaan we gebruik maken van de trackers die je net hebt toegewezen. Doe al je trackers aan, je kunt zien welke trackers welke zijn in de figuur rechts.
onboarding-automatic_mounting-put_trackers_on-next = Ik heb al mijn trackers aan

## Tracker manual proportions setupa

onboarding-manual_proportions-back = Ga terug naar de reset tutorial
onboarding-manual_proportions-title = Handmatige lichaamsverhoudingen
onboarding-manual_proportions-fine_tuning_button = Automatisch afstemmen van verhoudingen
onboarding-manual_proportions-fine_tuning_button-disabled-tooltip = Sluit een VR-headset aan om automatische fijnafstelling te gebruiken
onboarding-manual_proportions-export = Export proporties
onboarding-manual_proportions-import = Importeer proporties
onboarding-manual_proportions-file_type = Lichaamsproporties bestand
onboarding-manual_proportions-normal_increment = Normale verhoging
onboarding-manual_proportions-precise_increment = Nauwkeurige verhoging
onboarding-manual_proportions-grouped_proportions = Gegroepeerde verhoudingen
onboarding-manual_proportions-all_proportions = Alle verhoudingen
onboarding-manual_proportions-estimated_height = Geschatte gebruikerslengte

## Tracker automatic proportions setup

onboarding-automatic_proportions-back = Ga terug naar de reset tutorial
onboarding-automatic_proportions-title = Meet je lichaam
onboarding-automatic_proportions-description = Om SlimeVR-trackers te laten werken, moeten we de lengte van je botten weten. Deze korte kalibratie meet het voor je.
onboarding-automatic_proportions-manual = Handmatige kalibratie
onboarding-automatic_proportions-prev_step = Vorige stap
onboarding-automatic_proportions-put_trackers_on-title = Doe je trackers aan
onboarding-automatic_proportions-put_trackers_on-description = Om je verhoudingen te kalibreren, gaan we gebruik maken van de trackers die je net hebt toegewezen. Doe al je trackers aan, je kunt zien welke trackers welke zijn in de figuur rechts.
onboarding-automatic_proportions-put_trackers_on-next = Ik heb al mijn trackers aan
onboarding-automatic_proportions-requirements-title = Vereisten
# Each line of text is a different list item
onboarding-automatic_proportions-requirements-descriptionv2 = Je hebt voldaan aan de minimale vereisten om je voeten te tracken (over het algemeen 5 trackers). Je hebt je trackers en headset aan en draagt ze. Je trackers en headset zijn verbonden met de SlimeVR server en werken naar behoren (zonder haperingen, loskoppelingen etc.). Je headset stuurt positiedata naar de SlimeVR server (dit vereist doorgaans dat SteamVR draait en verbonden is met SlimeVR via de SlimeVR SteamVR-driver). De tracking werkt en registreert je bewegingen nauwkeurig (je hebt bijvoorbeeld een volledige reset uitgevoerd en de trackers bewegen in de juiste richting bij schoppen, bukken, zitten etc.).
onboarding-automatic_proportions-requirements-next = Ik heb de vereisten gelezen
onboarding-automatic_proportions-check_height-title-v3 = Meet de hoogte van uw headset
onboarding-automatic_proportions-check_height-description-v2 = De hoogte van uw headset (HMD) moet iets minder zijn dan uw volledige lengte, aangezien headsets uw ooghoogte meten. Deze meting wordt gebruikt als basis voor uw lichaamsverhoudingen.
# All the text is in bold!
onboarding-automatic_proportions-check_height-calculation_warning-v3 = Begin met meten terwijl je <u>rechtop</u> staat om je lengte te meten. Let erop dat je je handen niet hoger dan je headset tilt, want dat kan de meting beïnvloeden!
onboarding-automatic_proportions-check_height-guardian_tip = Als je een losse VR-bril gebruikt, zorg er dan voor dat je guardian/veilige zone is ingeschakeld zodat je lengte correct is gekalibreerd!
# Context is that the height is unknown
onboarding-automatic_proportions-check_height-unknown = Onbekend
# Shows an element below it
onboarding-automatic_proportions-check_height-hmd_height2 = De hoogte van uw headset is:
onboarding-automatic_proportions-check_height-measure-start = Begin met meten
onboarding-automatic_proportions-check_height-measure-stop = Stoppen met meten
onboarding-automatic_proportions-check_height-measure-reset = Probeer opnieuw te meten
onboarding-automatic_proportions-check_height-next_step = Ze zijn goed
onboarding-automatic_proportions-check_floor_height-title = Meet uw vloerhoogte (optioneel)
onboarding-automatic_proportions-check_floor_height-description = In sommige gevallen wordt uw vloerhoogte mogelijk niet correct ingesteld door uw headset, waardoor de hoogte van de headset hoger wordt gemeten dan zou moeten. U kunt de "hoogte" van uw vloer meten om de hoogte van uw headset te corrigeren.
# All the text is in bold!
onboarding-automatic_proportions-check_floor_height-calculation_warning-v2 = Begin met meten en zet een controller op je vloer om de hoogte te meten. Als je zeker weet dat je vloerhoogte klopt, kun je deze stap overslaan.
# Shows an element below it
onboarding-automatic_proportions-check_floor_height-floor_height = Uw vloerhoogte is:
onboarding-automatic_proportions-check_floor_height-full_height = Uw geschatte volledige lengte is:
onboarding-automatic_proportions-check_floor_height-measure-start = Begin met meten
onboarding-automatic_proportions-check_floor_height-measure-stop = Stoppen met meten
onboarding-automatic_proportions-check_floor_height-measure-reset = Probeer opnieuw te meten
onboarding-automatic_proportions-check_floor_height-skip_step = Sla deze stap over en sla op.
onboarding-automatic_proportions-check_floor_height-next_step = Gebruik vloerhoogte en bespaar
onboarding-automatic_proportions-start_recording-title = Zorg dat je klaar bent om te bewegen
onboarding-automatic_proportions-start_recording-description = We gaan nu enkele specifieke houdingen en bewegingen opnemen. Deze worden in het volgende scherm geprompt. Zorg dat je klaar bent om te beginnen als de knop wordt ingedrukt!
onboarding-automatic_proportions-start_recording-next = Start opname
onboarding-automatic_proportions-recording-title = REC
onboarding-automatic_proportions-recording-description-p0 = Opname bezig...
onboarding-automatic_proportions-recording-description-p1 = Voer de onderstaande bewegingen uit:
# Each line of text is a different list item
onboarding-automatic_proportions-recording-steps =
    Sta rechtop, rol je hoofd in een cirkel.
    Buig je rug naar voren en hurk. Kijk tijdens het hurken naar links en dan naar rechts.
    Draai je bovenlichaam naar links (tegen de klok in) en reik dan naar beneden naar de grond.
    Draai je bovenlichaam naar rechts (met de klok mee) en reik dan naar beneden naar de grond.
    Rol je heupen in een cirkelvormige beweging alsof je een hoelahoep gebruikt.
    Als er nog tijd over is voor de opname, kunt u deze stappen herhalen totdat deze is voltooid.
onboarding-automatic_proportions-recording-processing = Resultaat verwerken
# $time (Number) - Seconds left for the automatic calibration recording to finish (max 20)
onboarding-automatic_proportions-recording-timer =
    { $time ->
        [one] 1 seconde resterend
       *[other] { $time } seconden resterend
    }
onboarding-automatic_proportions-verify_results-title = Resultaten controleren
onboarding-automatic_proportions-verify_results-description = Controleer de resultaten hieronder, zien ze er correct uit?
onboarding-automatic_proportions-verify_results-results = Opnameresultaten
onboarding-automatic_proportions-verify_results-processing = Resultaat verwerken
onboarding-automatic_proportions-verify_results-redo = Opname opnieuw doen
onboarding-automatic_proportions-verify_results-confirm = Ze zijn correct
onboarding-automatic_proportions-done-title = Lichaam gemeten en opgeslagen.
onboarding-automatic_proportions-done-description = Je kalibratie voor lichaamsverhoudingen is voltooid!
onboarding-automatic_proportions-error_modal-v2 =
    <b>Waarschuwing:</b> Er is een fout opgetreden bij het schatten van de verhoudingen!
    Dit is waarschijnlijk een probleem met de montagekalibratie. Zorg ervoor dat je tracking goed werkt voordat je het opnieuw probeert.
     <docs>Bekijk de documentatie</docs> of word lid van onze <discord>Discord</discord> voor hulp ^_^
onboarding-automatic_proportions-error_modal-confirm = Begrepen!
onboarding-automatic_proportions-smol_warning =
    Uw ingestelde lengte van { $height } is lager dan de toegestane minimumlengte van { $minHeight }.
    <b>Voer de metingen opnieuw uit en controleer of ze correct zijn.</b>
onboarding-automatic_proportions-smol_warning-cancel = Ga terug

## Tracker scaled proportions setup

onboarding-scaled_proportions-title = Geschaalde proporties
onboarding-scaled_proportions-description =
    Voor een correcte werking van de SlimeVR-trackers hebben we de lengte van uw botten nodig.
    We gebruiken hiervoor een gemiddelde lichaamsverhouding, geschaald op basis van uw lengte.
onboarding-scaled_proportions-manual_height-title = Configureer uw lengte
onboarding-scaled_proportions-manual_height-description-v2 = Deze lengte wordt gebruikt als basis voor je lichaamsverhoudingen.
onboarding-scaled_proportions-manual_height-missing_steamvr =
    SteamVR is momenteel niet verbonden met SlimeVR, dus metingen kunnen niet worden gebaseerd op je headset.
    <b>Ga verder op eigen risico of raadpleeg de documentatie!</b>
onboarding-scaled_proportions-manual_height-height-v2 = Uw volledige lengte is
onboarding-scaled_proportions-manual_height-estimated_height = De geschatte hoogte van uw headset is:
onboarding-scaled_proportions-manual_height-next_step = Opslaan en doorgaan
onboarding-scaled_proportions-manual_height-warning =
    Je gebruikt momenteel de handmatige manier om geschaalde verhoudingen in te stellen!
    <b>Deze modus wordt alleen aanbevolen als je geen HMD met SlimeVR gebruikt.</b>
    
    Om de automatische geschaalde verhoudingen te kunnen gebruiken, doe het volgende:
onboarding-scaled_proportions-manual_height-warning-no_hmd = Sluit een VR-headset aan
onboarding-scaled_proportions-manual_height-warning-no_controllers = Zorg ervoor dat je controllers zijn verbonden en correct aan je handen zijn toegewezen

## Tracker scaled proportions reset

onboarding-scaled_proportions-reset_proportion-title = Reset je lichaamsverhoudingen
onboarding-scaled_proportions-reset_proportion-description = Om je lichaamsverhoudingen op basis van je lengte in te stellen, moet je nu al je verhoudingen resetten. Dit zal alle verhoudingen die je hebt ingesteld wissen en een basisconfiguratie bieden.
onboarding-scaled_proportions-done-title = Lichaamsverhoudingen ingesteld
onboarding-scaled_proportions-done-description = Je lichaamsverhoudingen zouden nu gebaseerd moeten zijn op je lengte

## Stay Aligned setup

onboarding-stay_aligned-title = Blijf in lijn
onboarding-stay_aligned-description = Stel Blijf in lijn in om je trackers in lijn te houden.
onboarding-stay_aligned-put_trackers_on-title = Doe je trackers aan
onboarding-stay_aligned-put_trackers_on-description = Om je ontspannen houdingen op te slaan, gebruiken we de trackers die je zojuist hebt toegewezen. Doe al je trackers om; je kunt zien welke welke zijn op de afbeelding rechts.
onboarding-stay_aligned-put_trackers_on-trackers_warning = Je hebt momenteel minder dan 5 trackers aangesloten en toegewezen! Dit is het minimale aantal trackers dat nodig is om Blijf in lijn goed te laten werken.
onboarding-stay_aligned-put_trackers_on-next = Ik heb al mijn trackers aan
onboarding-stay_aligned-verify_mounting-title = Controleer je montage
onboarding-stay_aligned-verify_mounting-step-0 = Blijf in lijn vereist een goede montage. Anders krijg je geen goede ervaring met Blijf in lijn.
onboarding-stay_aligned-verify_mounting-step-1 = 1. Beweeg terwijl je staat.
onboarding-stay_aligned-verify_mounting-step-2 = 2. Ga zitten en beweeg je benen en voeten.
onboarding-stay_aligned-verify_mounting-step-3 = 3. Als je trackers niet op de juiste plek zitten, druk dan op "Montagekalibratie opnieuw uitvoeren".
onboarding-stay_aligned-verify_mounting-redo_mounting = Montagekalibratie opnieuw uitvoeren
onboarding-stay_aligned-preparation-title = Voorbereiding
onboarding-stay_aligned-preparation-tip = Zorg dat je rechtop staat. Blijf recht vooruit kijken en houd je armen langs je zij.
onboarding-stay_aligned-relaxed_poses-standing-title = Ontspannen Staande Houding
onboarding-stay_aligned-relaxed_poses-standing-step-0 = 1. Ga in een comfortabele houding staan. Ontspan!
onboarding-stay_aligned-relaxed_poses-standing-step-1-v2 = 2. Druk op de knop "Houding opslaan".
onboarding-stay_aligned-relaxed_poses-sitting-title = Ontspannen Zittende Houding op Stoel
onboarding-stay_aligned-relaxed_poses-sitting-step-0 = 1. Ga comfortabel zitten. Ontspan!
onboarding-stay_aligned-relaxed_poses-sitting-step-1-v2 = 2. Druk op de knop "Houding opslaan".
onboarding-stay_aligned-relaxed_poses-flat-title = Ontspannen Zittende Houding op de Vloer
onboarding-stay_aligned-relaxed_poses-flat-step-0 = 1. Ga op de vloer zitten met je benen voor je uit. Ontspan!
onboarding-stay_aligned-relaxed_poses-flat-step-1-v2 = 2. Druk op de knop "Houding opslaan".
onboarding-stay_aligned-relaxed_poses-skip_step = Overslaan
onboarding-stay_aligned-done-title = Blijf in lijn ingeschakeld!
onboarding-stay_aligned-done-description = Je Blijf in lijn-instelling is voltooid!
onboarding-stay_aligned-done-description-2 = De setup is voltooid! Je kunt het proces opnieuw starten als je de houdingen opnieuw wilt kalibreren.
onboarding-stay_aligned-previous_step = Vorige
onboarding-stay_aligned-next_step = Volgende
onboarding-stay_aligned-restart = Herstarten
onboarding-stay_aligned-done = Klaar

## Home

home-no_trackers = Geen trackers gedetecteerd of toegewezen

## Trackers Still On notification

trackers_still_on-modal-title = Trackers staan nog steeds aan
trackers_still_on-modal-description =
    Een of meer trackers staan nog aan.
    Wil je SlimeVR toch afsluiten?
trackers_still_on-modal-confirm = SlimeVR afsluiten
trackers_still_on-modal-cancel = Wacht even...

## Status system

status_system-StatusTrackerReset = Het wordt aanbevolen om een volledige reset uit te voeren omdat een of meer trackers niet zijn aangepast.
status_system-StatusSteamVRDisconnected =
    { $type ->
        [steamvr_feeder] Momenteel niet verbonden naar de SlimeVR Feeder App.
       *[other] Momenteel niet verbonden naar SteamVR via de SlimeVR driver.
    }
status_system-StatusTrackerError = De { $trackerName } tracker heeft een error.
status_system-StatusUnassignedHMD = De VR-headset moet worden toegewezen als hoofdtracker.
status_system-StatusPublicNetwork =
    { $count ->
        [one] Je netwerkprofiel staat momenteel ingesteld op Openbaar ({ $adapters }). Dit wordt niet aanbevolen voor een goede werking van SlimeVR. <PublicFixLink>Hier lees je hoe je het kunt oplossen.</PublicFixLink>
       *[other] Sommige van je netwerkadapters staan ingesteld op openbaar: { $adapters }. Dit wordt niet aanbevolen voor een goede werking van SlimeVR. <PublicFixLink>Hier lees je hoe je dit kunt oplossen.</PublicFixLink>
    }

## Firmware tool globals

firmware_tool-next_step = Volgende stap
firmware_tool-previous_step = Vorige stap
firmware_tool-ok = Ziet er goed uit
firmware_tool-retry = Opnieuw
firmware_tool-loading = Laden...

## Firmware tool Steps

firmware_tool = DIY firmware-tool
firmware_tool-description = Hiermee kan je uw DIY-trackers configureren en flashen
firmware_tool-not_available = Oeps, de firmwaretool is momenteel niet beschikbaar. Kom later terug!
firmware_tool-not_compatible = De firmwaretool is niet compatibel met deze versie van de server. Gelieve te updaten!
firmware_tool-board_step = Selecteer je bord
firmware_tool-board_step-description = Selecteer een van de onderstaande borden.
firmware_tool-board_pins_step = Controleer de pinnen
firmware_tool-board_pins_step-description =
    Controleer of de geselecteerde pinnen correct zijn.
    Als je de SlimeVR-documentatie hebt gevolgd, zouden de standaardwaarden correct moeten zijn
firmware_tool-board_pins_step-enable_led = LED inschakelen
firmware_tool-board_pins_step-led_pin =
    .label = LED-pin
    .placeholder = Voer het adres van de LED-pin in
firmware_tool-board_pins_step-battery_type = Selecteer het batterijtype
firmware_tool-board_pins_step-battery_type-BAT_EXTERNAL = Externe batterij
firmware_tool-board_pins_step-battery_type-BAT_INTERNAL = Interne batterij
firmware_tool-board_pins_step-battery_type-BAT_INTERNAL_MCP3021 = Interne MCP3021
firmware_tool-board_pins_step-battery_type-BAT_MCP3021 = MCP3021
firmware_tool-board_pins_step-battery_sensor_pin =
    .label = Batterij sensor Pin
    .placeholder = Voer het pin-adres van de batterij sensor in
firmware_tool-board_pins_step-battery_resistor =
    .label = Batterij Weerstand (Ohm)
    .placeholder = Voer de waarde van de batterijweerstand in
firmware_tool-board_pins_step-battery_shield_resistor-0 =
    .label = Batterij Shield R1 (Ohm)
    .placeholder = Voer de waarde in van Battery Shield R1
firmware_tool-board_pins_step-battery_shield_resistor-1 =
    .label = Batterij Shield R2 (Ohm)
    .placeholder = Voer de waarde in van Battery Shield R2
firmware_tool-add_imus_step = Declareer uw IMU's
firmware_tool-add_imus_step-description =
    Voeg de IMU's toe die je tracker heeft
    Als je de SlimeVR-documentatie hebt gevolgd, zouden de standaardwaarden correct moeten zijn.
firmware_tool-add_imus_step-imu_type-label = IMU-type
firmware_tool-add_imus_step-imu_type-placeholder = Selecteer het type IMU
firmware_tool-add_imus_step-imu_rotation =
    .label = IMU-rotatie (graden)
    .placeholder = Rotatie van de IMU
firmware_tool-add_imus_step-scl_pin =
    .label = SCL-pin
    .placeholder = Pin-adres van SCL
firmware_tool-add_imus_step-sda_pin =
    .label = SDA-pin
    .placeholder = Pin-adres van SDA
firmware_tool-add_imus_step-int_pin =
    .label = INT-pin
    .placeholder = Pin-adres van INT
firmware_tool-add_imus_step-optional_tracker =
    .label = Optionele tracker
firmware_tool-add_imus_step-show_less = Toon minder
firmware_tool-add_imus_step-show_more = Toon meer
firmware_tool-add_imus_step-add_more = Voeg meer IMU's toe
firmware_tool-select_firmware_step = Selecteer de firmwareversie
firmware_tool-select_firmware_step-description = Kies de versie van de firmware die je wilt gebruiken
firmware_tool-select_firmware_step-show-third-party =
    .label = Firmware van derden weergeven
firmware_tool-flash_method_step = Flashing methode
firmware_tool-flash_method_step-description = Kies de flashingsmethode die je wilt gebruiken
firmware_tool-flash_method_step-ota =
    .label = OTA
    .description = Gebruik de draadloze methode. Je tracker zal de Wi-Fi gebruiken om de firmware bij te werken. Werkt alleen op reeds geconfigureerde trackers.
firmware_tool-flash_method_step-serial =
    .label = Serial
    .description = Gebruik een USB-kabel om je tracker bij te werken.
firmware_tool-flashbtn_step = Druk op de bootknop
firmware_tool-flashbtn_step-description = Voordat u naar de volgende stap gaat, zijn er een paar dingen die u moet doen.
firmware_tool-flashbtn_step-board_SLIMEVR = Zet de tracker uit, verwijder de behuizing (indien aanwezig), verbind een USB-kabel met deze computer en voer vervolgens een van de volgende stappen uit, afhankelijk van de revisie van uw SlimeVR-board:
firmware_tool-flashbtn_step-board_SLIMEVR-r11 = Zet de tracker aan terwijl u het tweede rechthoekige FLASH-pad vanaf de rand aan de bovenkant van het board kortsluit, en het metalen schild van de microcontroller.
firmware_tool-flashbtn_step-board_SLIMEVR-r12 = Zet de tracker aan terwijl u het ronde FLASH-pad aan de bovenkant van het board kortsluit, en het metalen schild van de microcontroller.
firmware_tool-flashbtn_step-board_SLIMEVR-r14 = Zet de tracker aan terwijl u de FLASH-knop aan de bovenkant van het board indrukt.
firmware_tool-flashbtn_step-board_OTHER =
    Voordat u gaat flashen, moet de tracker waarschijnlijk in de bootloader-modus worden gezet.
    Meestal betekent dit het indrukken van de bootknop op het board voordat het flashproces begint.
    Als het flashproces time-out bij het begin van het flashen, betekent dit waarschijnlijk dat de tracker niet in de bootloader-modus stond.
    Raadpleeg de flitsinstructies van uw board om te weten hoe u de bootloader-modus inschakelt.
firmware_tool-flash_method_ota-devices = Gedetecteerde OTA-apparaten:
firmware_tool-flash_method_ota-no_devices = Er zijn geen boards die via OTA bijgewerkt kunnen worden, zorg ervoor dat u het juiste boardtype heeft geselecteerd.
firmware_tool-flash_method_serial-wifi = Wi-Fi-gegevens:
firmware_tool-flash_method_serial-devices-label = Gedetecteerde serial apparaten:
firmware_tool-flash_method_serial-devices-placeholder = Selecteer een serieel apparaat
firmware_tool-flash_method_serial-no_devices = Er zijn geen compatibele seriële apparaten gedetecteerd, zorg ervoor dat de tracker is aangesloten.
firmware_tool-build_step = Aan het bouwen
firmware_tool-build_step-description = De firmware wordt gebouwd, even geduld a.u.b.
firmware_tool-flashing_step = Firmware aan het uploaden
firmware_tool-flashing_step-description = Je trackers worden geflashed, volg de instructies op het scherm
firmware_tool-flashing_step-warning-v2 = Koppel de tracker niet los en zet hem niet uit tijdens het uploadproces, tenzij dit wordt aangegeven. Dit kan je apparaat onbruikbaar maken.
firmware_tool-flashing_step-flash_more = Flash meer trackers
firmware_tool-flashing_step-exit = Sluit

## firmware tool build status

firmware_tool-build-CREATING_BUILD_FOLDER = De buildmap maken
firmware_tool-build-DOWNLOADING_FIRMWARE = Firmware wordt gedownload
firmware_tool-build-EXTRACTING_FIRMWARE = Firmware wordt uitgepakt
firmware_tool-build-SETTING_UP_DEFINES = Configureren van de definities
firmware_tool-build-BUILDING = Firmware wordt gebouwd
firmware_tool-build-SAVING = De build opslaan
firmware_tool-build-DONE = Build voltooid
firmware_tool-build-ERROR = Kan de firmware niet bouwen

## Firmware update status

firmware_update-status-DOWNLOADING = Firmware wordt gedownload
firmware_update-status-NEED_MANUAL_REBOOT-v2 = Zet je tracker uit en daarna weer aan.
firmware_update-status-AUTHENTICATING = Authenticatie met de mcu
firmware_update-status-UPLOADING = Firmware wordt geüpload
firmware_update-status-SYNCING_WITH_MCU = Synchroniseren met de mcu
firmware_update-status-REBOOTING = De update toepassen
firmware_update-status-PROVISIONING = Wi-Fi-inloggegevens instellen
firmware_update-status-DONE = Update voltooid!
firmware_update-status-ERROR_DEVICE_NOT_FOUND = Kan het apparaat niet vinden
firmware_update-status-ERROR_TIMEOUT = Er is een time-out opgetreden voor het updateproces
firmware_update-status-ERROR_DOWNLOAD_FAILED = Kan de firmware niet downloaden
firmware_update-status-ERROR_AUTHENTICATION_FAILED = Kan niet verifiëren met de mcu
firmware_update-status-ERROR_UPLOAD_FAILED = Kan de firmware niet uploaden
firmware_update-status-ERROR_PROVISIONING_FAILED = Kan de Wi-Fi-inloggegevens niet instellen
firmware_update-status-ERROR_UNSUPPORTED_METHOD = De updatemethode wordt niet ondersteund
firmware_update-status-ERROR_UNKNOWN = Onbekende fout

## Dedicated Firmware Update Page

firmware_update-title = Firmware-update
firmware_update-devices = Beschikbare apparaten
firmware_update-devices-description = Selecteer de trackers die u wilt updaten naar de nieuwste versie van SlimeVR-firmware
firmware_update-no_devices = Zorg er alsjeblieft voor dat de trackers die je wilt updaten AAN staan en verbonden zijn met de Wi-Fi!
firmware_update-changelog-title = Bijwerken naar { $version }
firmware_update-looking_for_devices = Op zoek naar apparaten om bij te werken...
firmware_update-retry = Opnieuw
firmware_update-update = Geselecteerde trackers bijwerken
firmware_update-exit = Sluit

## Tray Menu

tray_menu-show = Weergeven
tray_menu-hide = Verbergen
tray_menu-quit = Beëindigen

## First exit modal

tray_or_exit_modal-title = Wat is de actie van de sluitknop?
# Multiline text
tray_or_exit_modal-description =
    Hiermee kun je kiezen wat er gebeurt als je op de sluitknop klikt: het programma afsluiten of minimaliseren naar het systeemvak.
    
    Deze instelling kun je later altijd nog wijzigen in de interface instellingen!
tray_or_exit_modal-radio-exit = Afsluiten bij sluiten
tray_or_exit_modal-radio-tray = Minimaliseren naar systeemvak
tray_or_exit_modal-submit = Opslaan
tray_or_exit_modal-cancel = Annuleren

## Unknown device modal

unknown_device-modal-title = Er is een nieuwe tracker gevonden!
unknown_device-modal-description = Er is een nieuwe tracker gevonden met MAC-adres <b>{ $deviceId }</b>. Wil je deze verbinden met SlimeVR?
unknown_device-modal-confirm = Tuurlijk!
unknown_device-modal-forget = Negeer het
# VRChat config warnings
vrc_config-page-title = VRChat-configuratie waarschuwingen
vrc_config-page-desc = Deze pagina toont de status van je VRChat-instellingen en welke instellingen niet compatibel zijn met SlimeVR. Het wordt sterk aanbevolen om eventuele waarschuwingen hier te verhelpen voor de beste gebruikservaring met SlimeVR.
vrc_config-page-help = Kan je de instellingen niet finden?
vrc_config-page-help-desc = Bekijk onze <a>documentatie over dit onderwerp!</a>
vrc_config-page-big_menu = Tracking & IK (Hoofdmenu)
vrc_config-page-big_menu-desc = IK-instellingen in het hoofdmenu
vrc_config-page-wrist_menu = Tracking & IK (Polsmenu)
vrc_config-page-wrist_menu-desc = IK-instellingen in het kleine polsmenu
vrc_config-on = Aan
vrc_config-off = Uit
vrc_config-invalid = Je VRChat-instellingen zijn verkeerd geconfigureerd!
vrc_config-show_more = Toon meer
vrc_config-setting_name = Naam van de VRChat-instelling
vrc_config-recommended_value = Aanbevolen waarde
vrc_config-current_value = Huidige waarde
vrc_config-mute = Waarschuwing dempen
vrc_config-mute-btn = Dempen
vrc_config-unmute-btn = Dempen opheffen
vrc_config-legacy_mode = Use Legacy IK Solving
vrc_config-disable_shoulder_tracking = Disable Shoulder Tracking
vrc_config-shoulder_width_compensation = Shoulder Width Compensation
vrc_config-spine_mode = FBT Spine Mode
vrc_config-tracker_model = FBT Tracker Model
vrc_config-avatar_measurement_type = Avatar Measurement
vrc_config-calibration_range = Calibration Range
vrc_config-calibration_visuals = Display Calibration Visuals
vrc_config-user_height = User Real Height
vrc_config-spine_mode-UNKNOWN = Onbekend
vrc_config-spine_mode-LOCK_BOTH = Lock Both
vrc_config-spine_mode-LOCK_HEAD = Lock Head
vrc_config-spine_mode-LOCK_HIP = Lock Hip
vrc_config-tracker_model-UNKNOWN = Onbekend
vrc_config-tracker_model-AXIS = Axis
vrc_config-tracker_model-BOX = Box
vrc_config-tracker_model-SPHERE = Sphere
vrc_config-tracker_model-SYSTEM = Systeem
vrc_config-avatar_measurement_type-UNKNOWN = Onbekend
vrc_config-avatar_measurement_type-HEIGHT = Height
vrc_config-avatar_measurement_type-ARM_SPAN = Arm Span

## Error collection consent modal

error_collection_modal-title = Kunnen we fouten verzamelen?
error_collection_modal-description_v2 =
    { settings-interface-behavior-error_tracking-description_v2 }
    
    U kunt deze instelling later wijzigen in de sectie Gedrag van de instellingenpagina.
error_collection_modal-confirm = Ik ben akkoord
error_collection_modal-cancel = Ik wil het niet
